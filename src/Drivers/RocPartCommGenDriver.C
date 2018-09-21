#include <RocPartCommGenDriver.H>
#include <meshBase.H>
#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkUnstructuredGrid.h>
#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkGenericCell.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>
#include <vtkCleanPolyData.h>
#include <vtkAppendFilter.h>
#include <vtkPointLocator.h>
#include <vtkPointSet.h>
#include <vtkXMLPolyDataWriter.h>

#include <sstream>
#include <cstddef>
#include <AuxiliaryFunctions.H>
#include <cgnsWriter.H>
#include <iostream>
#include <fstream>
#include <typeinfo>

RocPartCommGenDriver::RocPartCommGenDriver(std::shared_ptr<meshBase> _mesh, 
                                           std::shared_ptr<meshBase> _remeshedSurf, 
                                           std::shared_ptr<meshBase> _volWithSol,
                                           std::shared_ptr<meshBase> _surfWithSol,
                                           std::shared_ptr<meshBase> _burnSurfWithSol,
                                           int numPartitions, const std::string& _base_t,
                                           bool _writeIntermediateFiles,
                                           double _searchTolerance, const std::string& _caseName)
{
  std::cout << "RocPartCommGenDriver created\n";
  this->mesh = _mesh;
  // load stitched surf mesh with patch info
  this->remeshedSurf = _remeshedSurf;
  this->volWithSol = _volWithSol;
  this->surfWithSol = _surfWithSol;
  this->burnSurfWithSol = _burnSurfWithSol;
  this->base_t = _base_t;
  this->caseName = _caseName;
  this->trimmed_base_t = std::to_string(std::stod(_base_t));
  this->trimmed_base_t.erase(this->trimmed_base_t.find_last_not_of('0')+1,std::string::npos);
  this->writeAllFiles = _writeIntermediateFiles;
  this->searchTolerance = _searchTolerance;
  this->execute(numPartitions);
}

RocPartCommGenDriver::RocPartCommGenDriver(std::string& volname, std::string& surfname, int numPartitions)
{
  std::cout << "RocPartCommGenDriver created\n";
  // load full volume mesh and create METIS partitions
  this->mesh = meshBase::CreateShared(volname);
  this->remeshedSurf = meshBase::CreateShared(surfname);
  this->base_t = "00.000000";
  this->trimmed_base_t = "0.0";
  this->writeAllFiles = false;
  this->execute(numPartitions);
}

void RocPartCommGenDriver::execute(int numPartitions)
{
  std::cout << "executing" << std::endl;
  remeshedSurf->setContBool(false);
  this->partitions = meshBase::partition(this->mesh.get(), numPartitions);
  this->AddGlobalCellIds(this->remeshedSurf);
  std::vector<std::string> paneDataAndGlobalCellIds = {"patchNo","bcflag","cnstr_type","GlobalCellIds"};
  // initialize storage for surf partitions 
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  this->surfacePartitions.resize(numPartitions);
  // initialize storage for virtual cells
  this->virtualCellsOfPartitions.resize(numPartitions);
  this->virtualCellsOfSurfPartitions.resize(numPartitions);
  // get all required local-global node maps and node/cell ids
  std::cout << "getting maps" << std::endl;
  this->getGlobalIdsAndMaps(numPartitions,true);
  // allocate storage for vol pconn vectors
  this->volPconns.resize(numPartitions);
  this->notGhostInPconn.resize(numPartitions);

  for (int i = 0; i < numPartitions; ++i)
  {
    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    appendFilter->AddInputData(
                    deleteInterPartitionSurface(
                      remeshedSurf, partitions[i]->extractSurface()));
    appendFilter->Update();
    // create surfs partitions casted from vtp to vtu
    unstructuredGrid = appendFilter->GetOutput();
    std::string basename("surfacePartition");
    basename += std::to_string(i);
    basename += ".vtu";
    // construct meshBase surface partition from vtkUnstructuredGrid
    this->surfacePartitions[i] =  meshBase::CreateShared(unstructuredGrid, basename);
    remeshedSurf->transfer(this->surfacePartitions[i].get(), 
                           "Consistent Interpolation", paneDataAndGlobalCellIds, 1);
    if (this->writeAllFiles) this->surfacePartitions[i]->write();
    // get ghost information for volume partitions
    this->getGhostInformation(i,true);
    // write pconn information for volume partition
    std::string type("01");
    notGhostInPconn[i] = this->writeSharedToPconn(i, type);
    this->writeSentToPconn(i,type,true);
    this->writeReceivedToPconn(i,type,true);
    this->writeSentToPconn(i,type,false);
    this->writeReceivedToPconn(i,type,false);

    this->comWriter(i+1);
    this->clearPconnVectors();

    for (int j = 0; j < numPartitions; ++j)
    {
      // get virtual cells of each volume partition (t4:virtual)
      this->getVirtualCells(i,j,true);
    }
    // write cgns for vol partition
    this->writeVolCgns("fluid", i,1,0);
    // write cell mapping file for entire mesh
    this->cmpWriter(0, this->mesh->getDataSet()->GetNumberOfCells());
    // write dimension file for entire mesh
    this->dimWriter(0, this->mesh, this->mesh);
  }
  // get global ids and maps for surface partitions
  this->getGlobalIdsAndMaps(numPartitions, false);
  for (int i = 0; i < numPartitions; ++i)
  {
    // get ghost information for each surface partition
    this->getGhostInformation(i,false);
    for (int j = 0; j < numPartitions; ++j)
    {
      // get virtual cells for each surface partition (t3:virtual)
      this->getVirtualCells(i,j,false); 
    }
  }
  // extract patches of each partition and get virtual cells from patches in other
  // partitions for each partiition
  this->extractPatches();
  // allocate shared patch nodes 4d map
  this->sharedPatchNodes.resize(numPartitions);
  // get shared nodes between patches both intra and inter partition
  this->getSharedPatchInformation();
  // allocate storage for surf pconn vectors
  this->surfPconns.resize(numPartitions);

  // Determining patch to zone numbers for surfaces
  // before writing surface Pconns
  std::map<int, int> surfZoneNumbersMap;
  for (int iProc = 0; iProc < numPartitions; ++iProc)
  {
    int zoneCounter = 2;
    surfZoneNumbersMap.clear();
    auto it = this->sharedPatchNodes[iProc].begin(); // patch
    while (it != sharedPatchNodes[iProc].end())
    {
      surfZoneNumbersMap[it->first] = zoneCounter;
      zoneCounter += 1;
      ++it;
    }
    surfZoneMap.push_back(surfZoneNumbersMap);
  }

  for (int i = 0; i < numPartitions; ++i)
  {
    this->writeSharedToPconn(i);
    this->writeSurfCgns("ifluid_b", i);
    this->writeSurfCgns("ifluid_ni", i);
    this->writeSurfCgns("burn", i);
    this->writeSurfCgns("iburn_all", i);
  }
  this->dimSurfWriter(0);

  this->txtWriter();

//=====================================================================================================
// !!!
// This is currently a really hacky way to accumulate the number of unique volume faces for the grid
// speed data entry. We have to accumulate the number of unique volume faces in first, then subtract
// off the surface faces, and then call the volume routine again to write out the actual total faces.
// (note the "gsExists" flag passed into writeVolCgns()). This works fine, but makes the driver run
// slower than it should.
// !!!
//=====================================================================================================

  remeshedSurf->setContBool(false);
  this->partitions = meshBase::partition(this->mesh.get(), numPartitions);
  this->AddGlobalCellIds(this->remeshedSurf);
  // initialize storage for surf partitions 
  this->surfacePartitions.resize(numPartitions);
  // initialize storage for virtual cells
  this->virtualCellsOfPartitions.resize(numPartitions);
  this->virtualCellsOfSurfPartitions.resize(numPartitions);
  // get all required local-global node maps and node/cell ids
  this->getGlobalIdsAndMaps(numPartitions,true);
  // allocate storage for vol pconn vectors
  this->volPconns.clear();
  this->volPconns.resize(numPartitions);
  this->notGhostInPconn.resize(numPartitions);
  for (int i = 0; i < numPartitions; ++i)
  {
    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    appendFilter->AddInputData(
                    deleteInterPartitionSurface(
                      remeshedSurf, partitions[i]->extractSurface()));
    appendFilter->Update();
    // create surfs partitions casted from vtp to vtu
    unstructuredGrid = appendFilter->GetOutput();
    std::string basename("surfacePartition");
    basename += std::to_string(i);
    basename += ".vtu";
    // construct meshBase surface partition from vtkUnstructuredGrid
    this->surfacePartitions[i] =  meshBase::CreateShared(unstructuredGrid, basename);
    remeshedSurf->transfer(this->surfacePartitions[i].get(), 
                           "Consistent Interpolation", paneDataAndGlobalCellIds, 1);
    if (this->writeAllFiles) this->surfacePartitions[i]->write();
    // get ghost information for volume partitions
    this->getGhostInformation(i,true);
    // write pconn information for volume partition
    std::string type("01");
    notGhostInPconn[i] = this->writeSharedToPconn(i, type);
    this->writeSentToPconn(i,type,true);
    this->writeReceivedToPconn(i,type,true);
    this->writeSentToPconn(i,type,false);
    this->writeReceivedToPconn(i,type,false);
    for (int j = 0; j < numPartitions; ++j)
    {
      // get virtual cells of each volume partition (t4:virtual)
      this->getVirtualCells(i,j,true);
    }
    // write cgns for vol partition
    this->writeVolCgns("fluid", i,1,1);
  }
//=====================================================================================================

  this->mapWriter();
}

void RocPartCommGenDriver::AddGlobalCellIds(std::shared_ptr<meshBase> _mesh)
{
  vtkSmartPointer<vtkDataArray> globalCellIds = vtkSmartPointer<vtkIdTypeArray>::New();
  globalCellIds->SetName("GlobalCellIds");
  globalCellIds->SetNumberOfComponents(1);
  globalCellIds->SetNumberOfTuples(_mesh->getNumberOfCells());
  for (int i = 0; i < _mesh->getNumberOfCells(); ++i)
  {
    globalCellIds->SetTuple1(i,i);
  }
  _mesh->getDataSet()->GetCellData()->AddArray(globalCellIds);  
}

void RocPartCommGenDriver::getGlobalIdsAndMaps(int numPartitions, bool vol)
{
  // clear vectors before use
  this->globToPartNodeMap.clear(); 
  this->globToPartCellMap.clear();
  this->partToGlobNodeMap.clear();
  this->partToGlobCellMap.clear();
  this->globalNodeIds.clear();
  this->globalCellIds.clear();
  // allocate
  this->globToPartNodeMap.resize(numPartitions);
  this->globToPartCellMap.resize(numPartitions);
  this->partToGlobNodeMap.resize(numPartitions);
  this->partToGlobCellMap.resize(numPartitions);
  this->globalNodeIds.resize(numPartitions);
  this->globalCellIds.resize(numPartitions);
  // get preloaded maps if volume
  if (vol)
  {
    for (int i = 0; i < numPartitions; ++i)
    {
      this->globToPartNodeMap[i] = partitions[i]->getGlobToPartNodeMap();
      this->globToPartCellMap[i] = partitions[i]->getGlobToPartCellMap();
      this->partToGlobNodeMap[i] = partitions[i]->getPartToGlobNodeMap();
      this->partToGlobCellMap[i] = partitions[i]->getPartToGlobCellMap();
      this->globalNodeIds[i] = getSortedKeys(this->globToPartNodeMap[i]);
      this->globalCellIds[i] = getSortedKeys(this->globToPartCellMap[i]);
    }
  }
  // compute maps and set if surface
  else
  {
    for (int i = 0; i < numPartitions; ++i)
    {
      vtkSmartPointer<vtkDataArray> globalNodeIds
        = surfacePartitions[i]->getDataSet()->GetPointData()->GetArray("GlobalNodeIds");
      for (int j = 0; j < surfacePartitions[i]->getNumberOfPoints(); ++j)
      { 
        int globNodeId = static_cast<int>(globalNodeIds->GetTuple1(j));
        this->globToPartNodeMap[i][globNodeId] = j;
        this->partToGlobNodeMap[i][j] = globNodeId;
      }
      vtkSmartPointer<vtkDataArray> globalCellIds
        = surfacePartitions[i]->getDataSet()->GetCellData()->GetArray("GlobalCellIds");
      for (int j = 0; j < surfacePartitions[i]->getNumberOfCells(); ++j)
      {
        int globCellId = static_cast<int>(globalCellIds->GetTuple1(j));
        this->globToPartCellMap[i][globCellId] = j;
        this->partToGlobCellMap[i][j] = globCellId;
      }
      this->globalNodeIds[i] = getSortedKeys(this->globToPartNodeMap[i]);
      this->globalCellIds[i] = getSortedKeys(this->globToPartCellMap[i]);
    }
  }  
}

int RocPartCommGenDriver::writeSharedToPconn(int proc, const std::string& type)
{
  auto sharedItr = this->sharedNodes[proc].begin();
  while ( sharedItr != this->sharedNodes[proc].end())
  {
    if (sharedItr->second.size() != 0)
    {
      // num blocks
      this->volPconns[proc].push_back(1);
      this->neighborProcsTmp.push_back(sharedItr->first);
      std::stringstream ss;
      ss << sharedItr->first + 1 << type;
      // zone number
      this->volPconns[proc].push_back(std::stoi(ss.str()));
      // number of ids
      this->volPconns[proc].push_back(sharedItr->second.size());
      // indices
      std::vector<int> tmpInd;
      for (int i = 0; i < sharedItr->second.size(); ++i)
      {
        this->volPconns[proc].push_back(sharedItr->second[i] + 1);
        tmpInd.push_back(sharedItr->second[i] + 1);
      }
      this->sharedNodesTmp.push_back(tmpInd);
    }
    ++sharedItr;
  }
  return volPconns[proc].size();
}

void RocPartCommGenDriver::writeSharedToPconn(int proc)
{
  auto it = this->sharedPatchNodes[proc].begin(); // patch
  while (it != sharedPatchNodes[proc].end())
  {
    auto it1 = it->second.begin();  // proc
    while (it1 != it->second.end())
    {
      auto it2 = it1->second.begin();  // patch
      while (it2 != it1->second.end())
      {
        // num blocks for this proc and current patch
        this->surfPconns[proc][it->first].push_back(1);
        // zone number
        std::stringstream ss;
        ss << it1->first + 1 << ((it->first+1) < 10 ? "0" : "") <<  this->surfZoneMap[it1->first][it2->first];
        this->surfPconns[proc][it->first].push_back(std::stoi(ss.str())); 
        // number of ids
        this->surfPconns[proc][it->first].push_back(it2->second.size());
        for (int k =0; k < it2->second.size(); ++k)
        {
          this->surfPconns[proc][it->first].push_back(it2->second[k]+1);
        } 
        ++it2;
      }
      ++it1;
    }
    ++it;
  }
}

void RocPartCommGenDriver::writeSentToPconn(int proc, const std::string& type, bool nodeOrCell)
{
  std::map<int, std::unordered_set<int>>::iterator sentItr;
  std::map<int, std::unordered_set<int>>::iterator end;
  if (nodeOrCell)
  {
    sentItr = this->sentNodes[proc].begin();
    end = this->sentNodes[proc].end();
  }
  else
  {
    sentItr = this->sentCells[proc].begin();
    end = this->sentCells[proc].end();
  }
  std::vector<int> tmpInd;
  while (sentItr != end)
  {
    if (sentItr->second.size() != 0)
    {
      // num blocks
      this->volPconns[proc].push_back(1);
      std::stringstream ss;
      ss << sentItr->first + 1 << type;
      // zone number
      this->volPconns[proc].push_back(std::stoi(ss.str()));
      // number of ids
      this->volPconns[proc].push_back(sentItr->second.size());
      // indices
      auto tmpIt = sentItr->second.begin();

      tmpInd.clear();
      while (tmpIt != sentItr->second.end())
      {
        this->volPconns[proc].push_back(*tmpIt + 1);
        tmpInd.push_back(*tmpIt + 1);
        ++tmpIt;
      }
      if (nodeOrCell)
      {
        this->sentNodesTmp.push_back(tmpInd);
      }
      else
      {
        this->sentCellsTmp.push_back(tmpInd);
      }
    } 
    ++sentItr;
  }
}

void RocPartCommGenDriver::writeReceivedToPconn(int proc, const std::string& type, bool nodeOrCell)
{
  std::map<int, std::unordered_set<int>>::iterator receivedItr;
  std::map<int, std::unordered_set<int>>::iterator end;
  int offset; 
  if (nodeOrCell)
  {
    receivedItr = this->receivedNodes[proc].begin();
    end = this->receivedNodes[proc].end();
    offset = partitions[proc]->getNumberOfPoints();
  }
  else
  {
    receivedItr = this->receivedCells[proc].begin();
    end = this->receivedCells[proc].end();
    offset = partitions[proc]->getNumberOfCells();
  }
  int numPrevReceived = 0;
  std::vector<int> tmpInd;
  while(receivedItr != end)
  {
    // num blocks  
    this->volPconns[proc].push_back(1);
    std::stringstream ss;
    ss << receivedItr->first + 1 << type;
    // zone number
    this->volPconns[proc].push_back(std::stoi(ss.str()));
    // number of ids; 
    int numReceived = receivedItr->second.size();
    this->volPconns[proc].push_back(numReceived);
    tmpInd.clear();
    for (int id = numPrevReceived; id < numPrevReceived+numReceived; ++id)
    {
      this->volPconns[proc].push_back(offset+id+1);
      tmpInd.push_back(offset+id+1);
    }
    if (nodeOrCell)
    {
      this->receivedNodesTmp.push_back(tmpInd);
    }
    else
    {
      this->receivedCellsTmp.push_back(tmpInd);
    }
    numPrevReceived += numReceived;
    ++receivedItr;
  }
}

void RocPartCommGenDriver::getVirtualCells(int me, int you, bool vol)
{
  if (vol && this->sharedNodes[me][you].size())
  {
    std::vector<int> virtuals(receivedCells[me][you].begin(), receivedCells[me][you].end());
    this->virtualCellsOfPartitions[me][you] = 
      meshBase::CreateShared(meshBase::extractSelectedCells(this->mesh.get(),virtuals));
    if (this->writeAllFiles)  
    {
      std::stringstream ss;
      ss << "virtual" << "Vol" << "Of" << me << "from" << you << ".vtu";
      this->virtualCellsOfPartitions[me][you]->write(ss.str());
    }
  }
  else if (!vol && this->sharedSurfNodes[me][you].size())
  {
    std::vector<int> virtuals(receivedSurfCells[me][you].begin(), receivedSurfCells[me][you].end());
    this->virtualCellsOfSurfPartitions[me][you] 
      = meshBase::CreateShared(meshBase::extractSelectedCells(this->remeshedSurf.get(),virtuals));
    if (this->writeAllFiles)  
    {
      std::stringstream ss;
      ss << "virtual" << "Surf" << "Of" << me << "from" << you << ".vtu";
      this->virtualCellsOfSurfPartitions[me][you]->write(ss.str());
    }
  }
}

void RocPartCommGenDriver::getGhostInformation(int me, bool volOrSurf)
{
  // will hold sent cell's indices
  vtkSmartPointer<vtkIdList> cellIdsList = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  // loop over all other procs to get nodes shared and nodes/cells sent
  for (int you = 0; you < partitions.size(); ++you)
  {
    if (you != me)
    {
      this->getGhostInformation(me,you, false, volOrSurf, cellIdsList, genCell);
      this->getGhostInformation(you,me, true, volOrSurf, cellIdsList, genCell);
    }
  }
}

void RocPartCommGenDriver::getGhostInformation(int me, int you, bool hasShared, bool vol,
                                               vtkSmartPointer<vtkIdList> cellIdsList,
                                               vtkSmartPointer<vtkGenericCell> genCell)
{
  // We need volume mesh to determine which ghost cells are included in the surface
  vtkSmartPointer<vtkIdList> volCellIdsList;
  meshBase* meMesh;
  meshBase* meVolMesh;
  if (vol)
  {
    meMesh = partitions[me].get();
  }
  else
  {
    meMesh = surfacePartitions[me].get();
    meVolMesh = partitions[me].get();
  }

  // ------ shared nodes between me and proc i
  // compute the intersection and assign to sharedNodes vec
  std::vector<int> tmpVec;
  std::set_intersection(globalNodeIds[me].begin(), globalNodeIds[me].end(),
                        globalNodeIds[you].begin(), globalNodeIds[you].end(),
                        std::back_inserter(tmpVec));

  if (!hasShared)
  {
    if (vol)
    {
      this->sharedNodes[me][you].resize(tmpVec.size());
      this->sharedNodes[you][me].resize(tmpVec.size());
    }
    else
    {
      this->sharedSurfNodes[me][you].resize(tmpVec.size());
      this->sharedSurfNodes[you][me].resize(tmpVec.size());
    }
  }
  // ------ fill shared nodes map and find cells and nodes me sends to proc i
  for (int j = 0; j < tmpVec.size(); ++j)
  {
    // get local idx of shared node in me
    int localPntId = globToPartNodeMap[me][tmpVec[j]];
    // add to sharedNodes maps for you and me
    if (!hasShared)
    {
      int procLocalPntId;
      if (vol)
      {
        this->sharedNodes[me][you][j] = localPntId;
        // get local idx of shared node in you
        procLocalPntId = globToPartNodeMap[you][tmpVec[j]];
        this->sharedNodes[you][me][j] = procLocalPntId; 
      }
      else
      {
        this->sharedSurfNodes[me][you][j] = localPntId;
        // get local idx of shared node in you
        procLocalPntId = globToPartNodeMap[you][tmpVec[j]];
        this->sharedSurfNodes[you][me][j] = procLocalPntId; 
      }
    }
    cellIdsList->Reset();
    // find cells using this shared node
    meMesh->getDataSet()->GetPointCells(localPntId, cellIdsList);
    // Nodes that are not shared among partitions
    std::vector<int> notSharedCellNodes;
    // Surface nodes that are not shared with the volume ghost cells's nodes
    std::vector<int> notSharedWithVolSurfNodes;
    // vtkIdList for checking if surface nodes are in volume mesh
    vtkSmartPointer<vtkIdList> result = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    // for each cell using shared node (these will be cells on the boundary)
    for (int k = 0; k < cellIdsList->GetNumberOfIds(); ++k)
    {
      // get the local cell idx
      int localCellId = cellIdsList->GetId(k);
      // get the cell in me for point extraction
      meMesh->getDataSet()->GetCell(localCellId, genCell);
      int numSharedInCell = 0;
      vtkSmartPointer<vtkPolyData> virtualNodesSet = vtkSmartPointer<vtkPolyData>::New();
      std::vector<int> virtualNodesIdList;
      if (!vol)
      {        
        // Construct a list of all volume ghost points. These will be compared against
        // the surface cells because ghost surface cells are a subset of ghost volume
        // cells.
        // TODO: This could be made faster.
        vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkPoints> virtualNodesList = vtkSmartPointer<vtkPoints>::New();
        // Construct a vtkDataSet of only the sentCells
        for (auto itr = sentCells[me][you].begin(); itr != sentCells[me][you].end(); ++itr)
        {
          meVolMesh->getDataSet()->GetCellPoints((*itr), cellPoints);
          for (int ipt = 0; ipt < cellPoints->GetNumberOfIds(); ipt++)
          {
            virtualNodesList->InsertNextPoint(meVolMesh->getDataSet()->GetPoint(cellPoints->GetId(ipt)));
            virtualNodesIdList.push_back(cellPoints->GetId(ipt));
          }
        }
        virtualNodesSet->SetPoints(virtualNodesList);
      }
      // Clear lists of nodes
      notSharedCellNodes.clear();
      notSharedWithVolSurfNodes.clear();
      // Keep track of virtual volume cell Ids that contain each potential virtual surface point Id
      std::vector<int> virtVolCellIds;
      virtVolCellIds.clear();
      for (int l = 0; l < genCell->GetNumberOfPoints(); ++l)
      {
        // get node idx of cell point
        int pntId = genCell->GetPointId(l);
        // if node is not found in shared nodes, it is sent 
        int globPntId = this->partToGlobNodeMap[me][pntId];
        if (std::find(tmpVec.begin(), tmpVec.end(), globPntId) == tmpVec.end())
        {
          if (vol)
          {
            notSharedCellNodes.push_back(pntId);
          }
          else
          {
            double* pntCoor = meMesh->getDataSet()->GetPoint(pntId);
            pointLocator->SetDataSet(virtualNodesSet);
            pointLocator->AutomaticOn();
            pointLocator->BuildLocator();
            // search tolerance
            double radius = 1e-8;
            result->Reset();
            pointLocator->FindPointsWithinRadius(radius, pntCoor, result);
            if (result->GetNumberOfIds())
            {
              // add idx to map of nodes me sends to you
              this->sentSurfNodes[me][you].insert(pntId);
              // add idx to map of nodes you recevies from proc me 
              this->receivedSurfNodes[you][me].insert(this->partToGlobNodeMap[me][pntId]);

              std::vector<int> virtVolCellIdsTmp;
              virtVolCellIdsTmp.clear();
              for (int i = 0; i < result->GetNumberOfIds(); i++)
              {
                vtkSmartPointer<vtkIdList> myList = vtkSmartPointer<vtkIdList>::New();
                meVolMesh->getDataSet()->GetPointCells(virtualNodesIdList[result->GetId(i)], myList);
                for (int i = 0; i < myList->GetNumberOfIds(); i++)
                {
                  virtVolCellIdsTmp.push_back(myList->GetId(i));
                }
              }
              // remove duplicates
              sort( virtVolCellIdsTmp.begin(), virtVolCellIdsTmp.end() );
              virtVolCellIdsTmp.erase( unique( virtVolCellIdsTmp.begin(), virtVolCellIdsTmp.end() ), virtVolCellIdsTmp.end() );
              for (auto itr = virtVolCellIdsTmp.begin(); itr != virtVolCellIdsTmp.end(); ++itr)
              {
                virtVolCellIds.push_back(*itr);
              }
            }
            else
            {
              notSharedWithVolSurfNodes.push_back(pntId);
            }
          }
        }
        else 
        {
          numSharedInCell += 1;

          if (!vol)
          {
            double* pntCoor = meMesh->getDataSet()->GetPoint(pntId);
            pointLocator->SetDataSet(virtualNodesSet);
            pointLocator->AutomaticOn();
            pointLocator->BuildLocator();
            // search tolerance
            double radius = 1e-8;
            result->Reset();
            pointLocator->FindPointsWithinRadius(radius, pntCoor, result);
            if (result->GetNumberOfIds())
            {
              // add idx to map of nodes me sends to you
              this->sentSurfNodes[me][you].insert(pntId);
              // add idx to map of nodes you recevies from proc me 
              this->receivedSurfNodes[you][me].insert(this->partToGlobNodeMap[me][pntId]);

              std::vector<int> virtVolCellIdsTmp;
              virtVolCellIdsTmp.clear();
              for (int i = 0; i < result->GetNumberOfIds(); i++)
              {
                vtkSmartPointer<vtkIdList> myList = vtkSmartPointer<vtkIdList>::New();
                meVolMesh->getDataSet()->GetPointCells(virtualNodesIdList[result->GetId(i)], myList);
                for (int i = 0; i < myList->GetNumberOfIds(); i++)
                {
                  virtVolCellIdsTmp.push_back(myList->GetId(i));
                }
              }
              // remove duplicates
              sort( virtVolCellIdsTmp.begin(), virtVolCellIdsTmp.end() );
              virtVolCellIdsTmp.erase( unique( virtVolCellIdsTmp.begin(), virtVolCellIdsTmp.end() ), virtVolCellIdsTmp.end() );
              for (auto itr = virtVolCellIdsTmp.begin(); itr != virtVolCellIdsTmp.end(); ++itr)
              {
                virtVolCellIds.push_back(*itr);
              }
            }
            else
            {
              notSharedWithVolSurfNodes.push_back(pntId);
            }
          }
        }
        if (vol && numSharedInCell >= 3)
        {
          // Add sent nodes only if the node belongs to a sent volume cell
          for (auto itr = notSharedCellNodes.begin(); itr != notSharedCellNodes.end(); itr++)
          {
            // add idx to map of nodes me sends to you
            this->sentNodes[me][you].insert(*itr);
            // add idx to map of nodes you recevies from proc me 
            this->receivedNodes[you][me].insert(this->partToGlobNodeMap[me][*itr]);
          }
        }
      }
      // sort list of accumulated virtual volume cell Ids that the potential virtual surface node
      // points belong to, and determine if all 3 (for tets) share a single cell
      int currId = -1;
      int accum = 0;
      int accumMax = 0;
      if (!vol)
      {
        std::sort(virtVolCellIds.begin(), virtVolCellIds.end());
        // Remove cells that aren't a part of the virtual sent cells
        std::vector<int> virtVolCellIdsTmp;
        for (auto itr = virtVolCellIds.begin(); itr != virtVolCellIds.end(); ++itr)
        {
          if (std::find(sentCells[me][you].begin(), sentCells[me][you].end(), *itr) != sentCells[me][you].end())
          {
            virtVolCellIdsTmp.push_back(*itr);
          }
        }
        virtVolCellIds = virtVolCellIdsTmp;

        for (auto itr = virtVolCellIds.begin(); itr != virtVolCellIds.end(); ++itr)
        {
          if (*itr != currId)
          {
            currId = *itr;
            accum = 0;
          }
          accum += 1;
          if (accum > accumMax)
          {
            accumMax = accum;
          }
        }
      }
      if (vol && numSharedInCell >= 3)
      {
        // add idx to sentCells to you map
        this->sentCells[me][you].insert(localCellId);
        // get the global cell index
        int globId = partToGlobCellMap[me][localCellId];
        // add idx to map of cells proc you recieves from proc me
        this->receivedCells[you][me].insert(globId);
      }
      else if (!vol && numSharedInCell >= 2 && notSharedWithVolSurfNodes.empty() && accumMax == 3)
      {
        int globId = partToGlobCellMap[me][localCellId];
        this->sentSurfCells[me][you].insert(localCellId);
        // add idx to map of cells proc you recieves from proc me
        this->receivedSurfCells[you][me].insert(globId);
      }
    }
  }
}

void RocPartCommGenDriver::getSharedPatchInformation()
{
  // get glob to part node maps for intersection calc
  std::map<int, std::map<int, std::map<int,int>>> patchProcGlobToPartNodeMaps;
  std::map<int, std::map<int, std::map<int,int>>> patchProcPartToGlobNodeMaps;
  std::map<int, std::map<int, std::vector<int>>> patchProcGlobNodeIdsMaps;

  for (int i = 0; i < this->sharedPatchNodes.size(); ++i)
  {
    auto it = this->patchesOfSurfacePartitions[i].begin();
    while (it != this->patchesOfSurfacePartitions[i].end())
    {
      // get global ids of patch it->first on proc i (only if patch mesh exists)
      vtkSmartPointer<vtkDataArray> vtkGlobPatchNodeIds
        = it->second->getDataSet()->GetPointData()->GetArray("GlobalNodeIds");
      std::map<int,int> patchGlobToPartNodeMap;
      std::map<int,int> patchPartToGlobNodeMap;
      for (int j = 0; j < vtkGlobPatchNodeIds->GetNumberOfTuples(); ++j)
      {
        int globNodeId = static_cast<int>(vtkGlobPatchNodeIds->GetTuple1(j));
        patchGlobToPartNodeMap[globNodeId] = j;
        patchPartToGlobNodeMap[j] = globNodeId;
      }
      patchProcGlobToPartNodeMaps[it->first][i] = patchGlobToPartNodeMap;
      patchProcPartToGlobNodeMaps[it->first][i] = patchPartToGlobNodeMap;
      patchProcGlobNodeIdsMaps[it->first][i] = getSortedKeys(patchGlobToPartNodeMap);
      ++it;
    }
  }
  for (int me = 0; me < this->surfacePartitions.size(); ++me)
  {
    auto it = this->patchesOfSurfacePartitions[me].begin();
    while (it != this->patchesOfSurfacePartitions[me].end())
    {   
      if (!(patchProcGlobNodeIdsMaps[it->first].find(me) == patchProcGlobNodeIdsMaps[it->first].end()))
      {  
        for (int you = 0; you < this->surfacePartitions.size(); ++you)
        {
          auto it1 = this->patchesOfSurfacePartitions[you].begin();
          while (it1 != this->patchesOfSurfacePartitions[you].end())
          {  
            if (!(patchProcGlobNodeIdsMaps[it1->first].find(you) == patchProcGlobNodeIdsMaps[it1->first].end()))
            {
              if (!(you == me && it->first == it1->first))
              {
                std::vector<int> tmpVec;
                std::set_intersection(patchProcGlobNodeIdsMaps[it->first][me].begin(), 
                                      patchProcGlobNodeIdsMaps[it->first][me].end(),
                                      patchProcGlobNodeIdsMaps[it1->first][you].begin(),
                                      patchProcGlobNodeIdsMaps[it1->first][you].end(),
                                      std::back_inserter(tmpVec));
                if (tmpVec.size())
                {
                  this->sharedPatchNodes[me][it->first][you][it1->first].resize(tmpVec.size());
                  this->sharedPatchNodes[you][it1->first][me][it->first].resize(tmpVec.size());
                  for (int i = 0; i < tmpVec.size(); ++i)
                  {
                    int localPntId = patchProcGlobToPartNodeMaps[it->first][me][tmpVec[i]];
                    this->sharedPatchNodes[me][it->first][you][it1->first][i] = localPntId;
                    localPntId = patchProcGlobToPartNodeMaps[it1->first][you][tmpVec[i]];
                    this->sharedPatchNodes[you][it1->first][me][it->first][i] = localPntId;
                  }
                }
              }
            }
            ++it1;
          }                    
        }
      }
      ++it;
    }
  }
}

vtkSmartPointer<vtkPolyData>
RocPartCommGenDriver::deleteInterPartitionSurface(std::shared_ptr<meshBase> fullSurf, 
                                              vtkSmartPointer<vtkDataSet> _partSurf)
{
  vtkSmartPointer<vtkPolyData> partSurf
    = vtkPolyData::SafeDownCast(_partSurf);
  // build cell locator for full surface
  vtkSmartPointer<vtkCellLocator> fullSurfCellLocator = fullSurf->buildLocator();
  // build upward links from points to cells
  partSurf->BuildLinks();
  // allocate generic cell for FindCell calls
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
  std::unordered_set<int> cellsToDelete;
  std::unordered_set<int> pointsToDelete;
  int numToDelete = 0;

  for (int i = 0; i < partSurf->GetNumberOfCells(); ++i)
  {
    // get center of cell from part surf
    double center[3] = {0.,0.,0.};
    partSurf->GetCell(i,genCell);
    for (int j = 0; j < genCell->GetNumberOfPoints(); ++j)
    {
      double point[3];
      genCell->GetPoints()->GetPoint(j,point);
      for (int k = 0; k < 3; ++k)
      {
        center[k] += point[k]/3.;
      }  
    }
    // params to find center in full surf with locator
    double closestPoint[3];
    vtkIdType id;
    int subid;
    double minDist2;
    // look for point in full surf cells
    fullSurfCellLocator->FindClosestPoint(center, closestPoint, genCell, id, subid, minDist2);  
    if (minDist2 > this->searchTolerance)
    {
      cellsToDelete.insert(i);
      partSurf->GetCellPoints(i, cellPointIds); 
      for (int l = 0; l < cellPointIds->GetNumberOfIds(); ++l)
      {
        int pntId = static_cast<int>(cellPointIds->GetId(l));
        pointsToDelete.insert(pntId);
      }
      cellPointIds->Reset();
      numToDelete += 1;
    }
  } 
  std::unordered_set<int>::iterator it = pointsToDelete.begin();
  while (it != pointsToDelete.end())
  {
    partSurf->DeletePoint(*it);
    ++it;
  }
  it = cellsToDelete.begin();
  while (it != cellsToDelete.end())
  {
    partSurf->DeleteCell(*it);
    ++it;
  }
  partSurf->RemoveDeletedCells();
  // force removal of cells, points and duplicates
  vtkSmartPointer<vtkCleanPolyData> cleanPolyData 
    = vtkSmartPointer<vtkCleanPolyData>::New();
  cleanPolyData->SetInputData(partSurf);
  cleanPolyData->Update();
  return cleanPolyData->GetOutput();
}

void RocPartCommGenDriver::extractPatches()
{
  this->patchesOfSurfacePartitions.resize(surfacePartitions.size());
  for (int i = 0; i < this->surfacePartitions.size(); ++i)
  {
    std::map<int, std::vector<int>> patchPartitionCellMap;
    vtkSmartPointer<vtkDataArray> patchNumbers
      = this->surfacePartitions[i]->getDataSet()->GetCellData()->GetArray("patchNo");
    for (int j = 0; j < patchNumbers->GetNumberOfTuples(); ++j)
    {
      int patchNo = static_cast<int>(patchNumbers->GetTuple1(j));
      patchPartitionCellMap[patchNo].push_back(j);
    }
    auto it = patchPartitionCellMap.begin();
    while (it != patchPartitionCellMap.end())
    {
      this->patchesOfSurfacePartitions[i][it->first] 
        = meshBase::CreateShared(
            meshBase::extractSelectedCells(this->surfacePartitions[i].get(), it->second));
      if (this->writeAllFiles)
      {
        std::stringstream ss;
        ss << "extractedPatch" << it->first << "OfProc" << i << ".vtu";
        this->patchesOfSurfacePartitions[i][it->first]->write(ss.str());
      }
      ++it;
    }
  }
  this->virtualCellsOfPatchesOfSurfacePartitions.resize(surfacePartitions.size());
  for (int i = 0; i < this->virtualCellsOfSurfPartitions.size(); ++i)
  {
    auto it = virtualCellsOfSurfPartitions[i].begin();
    while (it != virtualCellsOfSurfPartitions[i].end())
    {
      std::map<int, std::vector<int>> virtualPatchPartitionCellMap;
      vtkSmartPointer<vtkDataArray> virtualPatchNumbers
        = it->second->getDataSet()->GetCellData()->GetArray("patchNo");
      for (int j = 0; j < virtualPatchNumbers->GetNumberOfTuples(); ++j)
      {
        int patchNo = static_cast<int>(virtualPatchNumbers->GetTuple1(j));
        virtualPatchPartitionCellMap[patchNo].push_back(j);
      }
      auto it1 = virtualPatchPartitionCellMap.begin();
      while (it1 != virtualPatchPartitionCellMap.end())
      {
        this->virtualCellsOfPatchesOfSurfacePartitions[i][it->first][it1->first]
          = meshBase::CreateShared(
              meshBase::extractSelectedCells(it->second.get(), it1->second));
        if (this->writeAllFiles)
        {
          std::stringstream ss;
          ss << "extractedVirtualPatch" << it1->first 
            << "Of" << i << "FromProc" << it->first  << ".vtu";
          this->virtualCellsOfPatchesOfSurfacePartitions[i][it->first][it1->first]->write(ss.str());
        }
        ++it1;
      }
      ++it;
    }   
  }
}

RocPartCommGenDriver* RocPartCommGenDriver::readJSON(json inputjson)
{
  std::string volname = inputjson["Remeshed Volume"].as<std::string>();
  std::string surfname = inputjson["Stitched Surface"].as<std::string>();
  int numPartitions = inputjson["Number Of Partitions"].as<int>();
  return new RocPartCommGenDriver(volname, surfname, numPartitions);
}

RocPartCommGenDriver::~RocPartCommGenDriver()
{
  std::cout << "RocPartCommGenDriver destroyed" << std::endl;
}

void RocPartCommGenDriver::writeSurfCgns(const std::string& prefix, int me)
{
  std::stringstream ss;
  ss << prefix << "_" << this->base_t << 
    (me < 1000 ? (me < 100 ? (me < 10 ? "_000" : "_00") : "_0") : "_") << me << ".cgns";
  std::string cgFname(ss.str());
  std::unique_ptr<cgnsWriter> writer
    = std::unique_ptr<cgnsWriter>(new cgnsWriter(cgFname,prefix,2,3));
  // initialize Rocstar units, dimensions, and magnitude writing
  writer->setFluidUnitsMap();
  writer->setFluidDimMap();
  writer->setFluidMagMap();
  writer->setiFluidUnitsMap();
  writer->setiFluidDimMap();
  writer->setiFluidMagMap();
  writer->setBurnUnitsMap();
  writer->setBurnDimMap();
  writer->setBurnMagMap();
  // set time stamp
  writer->setTimestamp(this->trimmed_base_t);
  // define elementary information
  writer->setUnits(Kilogram, Meter, Second, Kelvin, Degree);
  // baseitrname, nTstep, timeVal
  writer->setBaseItrData("TimeIterValues", 1, std::stod(this->trimmed_base_t));
  // write elementary grid information to file
  writer->writeGridToFile();
  std::string grdpntr("Grid");
  grdpntr += this->trimmed_base_t;
  std::string flowsolpntr("NodeData");
  flowsolpntr += this->trimmed_base_t;
  ss.str("");
  ss.clear();
  if (!(prefix == "burn" || prefix == "iburn_all"))
  {
    // set Win data
    std::string winName("WinData");
    winName += this->trimmed_base_t;
    writer->setIntData(winName, 1);
    writer->writeWinToFile();
  }

  // Always one volumetric "fluid" zone (1)
  // Currently assumes exactly 1 burning patch
  // TODO: Support more than 1 burning patch
  // TODO: add support for "nb" files (currently only "b" and "ni")
  int zone = 2;
  if (this->patchesOfSurfacePartitions[me].find(1) 
        != this->patchesOfSurfacePartitions[me].end() && prefix == "ifluid_ni")
  {
    zone += 1; // zone starts at 3 in case there is a burning patch on partition 
  }

  // vtkIdList for mapping surface node indices to volume indices
  vtkSmartPointer<vtkIdList> result = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
  vtkSmartPointer<vtkPoints> myPoints = vtkSmartPointer<vtkPoints>::New();
  for (vtkIdType i = 0; i < procVolMesh[me]->getDataSet()->GetNumberOfPoints(); i++)
  {
    myPoints->InsertNextPoint(procVolMesh[me]->getDataSet()->GetPoint(i));
  }
  vtkSmartPointer<vtkPolyData> myPointsData = vtkSmartPointer<vtkPolyData>::New();
  myPointsData->SetPoints(myPoints);
  pointLocator->SetDataSet(myPointsData);
  pointLocator->AutomaticOn();
  pointLocator->BuildLocator();

  // search tolerance for VTK searches
  double radius = 1e-8;

  // begin loop over patches of this partition
  auto it = this->patchesOfSurfacePartitions[me].begin();
  int patchesFound = 0;
  while (it != this->patchesOfSurfacePartitions[me].end())
  {
    if ((prefix == "ifluid_ni" && it->first != 1) || (prefix == "ifluid_b" && it->first == 1) ||
          (prefix == "burn" && it->first == 1) || (prefix == "iburn_all" && it->first == 1)) // ni files have zones > 1 and 1 is a burning patch number
    {
      patchesFound += 1;
      int numVirtualCells = 0;
      // setting zone iter data
      writer->setZoneItrData("ZoneIterativeData", grdpntr, flowsolpntr);
      ss.str("");
      ss.clear();
      //ss << 0 << me+1 << (it->first < 10 ? "0" : "") <<  zone;
      ss << ((me+1) >= 10 ? "" : "0") << me+1 << (it->first < 10 ? "0" : "") << zone;
      zone += 1;    
      // vector to hold real patch mesh and all virtuals from same patch on other procs
      std::vector<std::shared_ptr<meshBase>> patchOfPartitionWithAllVirtualCells;
      // vector to hold real patch mesh and all virtuals from same patch on other procs
      std::vector<std::shared_ptr<meshBase>> patchOfPartition;
      // first push back real cells
      patchOfPartitionWithAllVirtualCells.push_back(it->second);
      patchOfPartition.push_back(it->second);
      // iterate over virtual cells of the partition's patch and push back
      auto virtItr = this->virtualCellsOfPatchesOfSurfacePartitions[me].begin();
      while (virtItr != this->virtualCellsOfPatchesOfSurfacePartitions[me].end())
      {
        auto virtItr1 = virtItr->second.begin();
        while (virtItr1 != virtItr->second.end())
        {
          // if the patch is the same as the one on the other proc
          if (it->first == virtItr1->first) 
          {
            if (virtItr1->second->getNumberOfCells())  
            {
              // push back virtual cells from same patch on other proc
              patchOfPartitionWithAllVirtualCells.push_back(virtItr1->second);
              numVirtualCells += virtItr1->second->getNumberOfCells(); 
            }
          }
          ++virtItr1;
        } 
        ++virtItr;
      } 

      if (numVirtualCells || !numVirtualCells)
      { 
        // merge real cells for this partition
        std::shared_ptr<meshBase> patchOfPartitionWithRealMesh
          = meshBase::stitchMB(patchOfPartition);
        // merge virtual cells into real mesh for this partition
        std::shared_ptr<meshBase> patchOfPartitionWithVirtualMesh
          = meshBase::stitchMB(patchOfPartitionWithAllVirtualCells);
        // set the zone
        writer->setZone(ss.str(), Unstructured);
        // all pconn stuff is shared nodes, so no ghost entities
        writer->setPconnGhostDescriptor(0);
        //writer->setPconnGhostDescriptor(numVirtualCells);
        // set num real vertices
        writer->setNVrtx(it->second->getNumberOfPoints());
        // set num real cells for this patch of this partition
        writer->setNCell(it->second->getNumberOfCells());
        // define and set coordinates
        std::vector<std::vector<double>> coords;

        if (!(prefix == "burn" || prefix == "iburn_all"))
        {
          coords = patchOfPartitionWithVirtualMesh->getVertCrds();
        }
        else
        {
          coords = patchOfPartitionWithVirtualMesh->getVertCrds();
        }

        // define real connectivities and 1-base indexing
        double myPoint[3];
        for (vtkIdType i = 0; i < it->second->getDataSet()->GetNumberOfPoints(); i++)
        {
            it->second->getDataSet()->GetPoint(i, myPoint);
        }
        std::vector<int> cgConnReal(it->second->getConnectivities());
        for (auto tmpit = cgConnReal.begin(); tmpit != cgConnReal.end(); ++tmpit)
        {
          *tmpit += 1;
        }
        // define virtual connectivies and 1-base indexing
        std::vector<int> cgConnVirtual(patchOfPartitionWithVirtualMesh->getConnectivities());
        cgConnVirtual.erase(cgConnVirtual.begin(), cgConnVirtual.begin()+cgConnReal.size());

        for (auto tmpit = cgConnVirtual.begin(); tmpit != cgConnVirtual.end(); ++tmpit)
        {
          *tmpit += 1;
        }
        this->swapTriOrdering(cgConnReal);
        this->swapTriOrdering(cgConnVirtual);

        // Construct surface-to-volume mapping of vertex IDs
        // Create vectors for t3g indices
        std::vector<int> cgConnRealGlobal1;
        std::vector<int> cgConnRealGlobal2;
        std::vector<int> cgConnRealGlobal3;
        std::vector<int> cgConnVirtualGlobal1;
        std::vector<int> cgConnVirtualGlobal2;
        std::vector<int> cgConnVirtualGlobal3;
        // Map real vertices
        vtkIdType foundId;
        double pntCoor[3];
        //Construct map from local to global vertices and vice versa
        std::map<int,int> loc2Glob;
        std::map<int,int> glob2Loc;
        // Get Global Real IDs
        for (auto itr = cgConnReal.begin(); itr != cgConnReal.end(); ++itr)
        {
          // Undo 1-base indexing
          pntCoor[0] = coords[0][*itr-1];
          pntCoor[1] = coords[1][*itr-1];
          pntCoor[2] = coords[2][*itr-1];
          foundId = pointLocator->FindClosestPoint(pntCoor);
          cgConnRealGlobal1.push_back(foundId+1); // redo 1-base indexing
          loc2Glob[*itr] = foundId+1;   // build up local-to-global map
          glob2Loc[foundId+1] = *itr;   // build up global-to-local map

          ++itr;
          pntCoor[0] = coords[0][*itr-1];
          pntCoor[1] = coords[1][*itr-1];
          pntCoor[2] = coords[2][*itr-1];
          foundId = pointLocator->FindClosestPoint(pntCoor);
          cgConnRealGlobal2.push_back(foundId+1);
          loc2Glob[*itr] = foundId+1;
          glob2Loc[foundId+1] = *itr;

          ++itr;
          pntCoor[0] = coords[0][*itr-1];
          pntCoor[1] = coords[1][*itr-1];
          pntCoor[2] = coords[2][*itr-1];
          foundId = pointLocator->FindClosestPoint(pntCoor);
          cgConnRealGlobal3.push_back(foundId+1);
          loc2Glob[*itr] = foundId+1;
          glob2Loc[foundId+1] = *itr;
        }
        // Get Global Virtual IDs
        if (numVirtualCells)
        {
          for (auto itr = cgConnVirtual.begin(); itr != cgConnVirtual.end(); ++itr)
          {
            pntCoor[0] = coords[0][*itr-1];
            pntCoor[1] = coords[1][*itr-1];
            pntCoor[2] = coords[2][*itr-1];
            foundId = pointLocator->FindClosestPoint(pntCoor);
            cgConnVirtualGlobal1.push_back(foundId+1);
            if (loc2Glob.find(*itr) == loc2Glob.end())
            {
              loc2Glob[*itr] = foundId+1;
              glob2Loc[foundId+1] = *itr;
            }
            double myPoint[3];
            pointLocator->GetDataSet()->GetPoint(foundId, myPoint);
            ++itr;

            pntCoor[0] = coords[0][*itr-1];
            pntCoor[1] = coords[1][*itr-1];
            pntCoor[2] = coords[2][*itr-1];
            foundId = pointLocator->FindClosestPoint(pntCoor);
            cgConnVirtualGlobal2.push_back(foundId+1);
            if (loc2Glob.find(*itr) == loc2Glob.end())
            {
              loc2Glob[*itr] = foundId+1;
              glob2Loc[foundId+1] = *itr;
            }
            pointLocator->GetDataSet()->GetPoint(foundId, myPoint);

            ++itr;
            pntCoor[0] = coords[0][*itr-1];
            pntCoor[1] = coords[1][*itr-1];
            pntCoor[2] = coords[2][*itr-1];
            foundId = pointLocator->FindClosestPoint(pntCoor);
            cgConnVirtualGlobal3.push_back(foundId+1);
            if (loc2Glob.find(*itr) == loc2Glob.end())
            {
              loc2Glob[*itr] = foundId+1;
              glob2Loc[foundId+1] = *itr;
            }
            pointLocator->GetDataSet()->GetPoint(foundId, myPoint);
          }
        }
        else
        {
          cgConnVirtualGlobal1.push_back(0);
          cgConnVirtualGlobal2.push_back(0);
          cgConnVirtualGlobal3.push_back(0);
        }

        // Maps between old and new local Ids
        std::map<int,int> old2New;
        std::map<int,int> new2Old;
        if (numVirtualCells || !numVirtualCells)
        {
          // Sort coordinates according to global index
          std::vector<std::vector<double>> coordsTmp;
          std::vector<double> coordsXTmp;
          std::vector<double> coordsYTmp;
          std::vector<double> coordsZTmp;
          int gId, lId;
          for (auto itr = glob2Loc.begin(); itr != glob2Loc.end(); ++itr)
          {
            gId = itr->first;
            lId = itr->second;
            coordsXTmp.push_back(coords[0][lId-1]);
            coordsYTmp.push_back(coords[1][lId-1]);
            coordsZTmp.push_back(coords[2][lId-1]);
          }
          coordsTmp.push_back(coordsXTmp);
          coordsTmp.push_back(coordsYTmp);
          coordsTmp.push_back(coordsZTmp);        

          // Sort local and global real connectivites
          std::vector<int> cgConnRealTmp;
          std::vector<int> cgConnRealGlobal1Tmp;
          std::vector<int> cgConnRealGlobal2Tmp;
          std::vector<int> cgConnRealGlobal3Tmp;
          int oldLid;
          int newLid;

          for (auto itr = cgConnReal.begin(); itr != cgConnReal.end(); ++itr)
          {
            oldLid = *itr;
            newLid = std::distance(glob2Loc.begin(), glob2Loc.find(loc2Glob[*itr])) + 1;
            old2New[oldLid] = newLid;
            new2Old[newLid] = oldLid;
            cgConnRealTmp.push_back(newLid);
            cgConnRealGlobal1Tmp.push_back(loc2Glob[oldLid]);
            ++itr;
  
            oldLid = *itr;
            newLid = std::distance(glob2Loc.begin(), glob2Loc.find(loc2Glob[*itr])) + 1;
            old2New[oldLid] = newLid;
            new2Old[newLid] = oldLid;
            cgConnRealTmp.push_back(newLid);
            cgConnRealGlobal2Tmp.push_back(loc2Glob[oldLid]);
            ++itr;
  
            oldLid = *itr;
            newLid = std::distance(glob2Loc.begin(), glob2Loc.find(loc2Glob[*itr])) + 1;
            old2New[oldLid] = newLid;
            new2Old[newLid] = oldLid;
            cgConnRealTmp.push_back(newLid);
            cgConnRealGlobal3Tmp.push_back(loc2Glob[oldLid]);
          }
  
          // Sort local and global virtual connectivites
          std::vector<int> cgConnVirtualTmp;
          std::vector<int> cgConnVirtualGlobal1Tmp;
          std::vector<int> cgConnVirtualGlobal2Tmp;
          std::vector<int> cgConnVirtualGlobal3Tmp;
          int oldGid;
          int newGid;

          if (numVirtualCells)
          {
            for (auto itr = cgConnVirtual.begin(); itr != cgConnVirtual.end(); ++itr)
            {
              oldLid = *itr;
              newLid = std::distance(glob2Loc.begin(), glob2Loc.find(loc2Glob[*itr])) + 1;
              old2New[oldLid] = newLid;
              new2Old[newLid] = oldLid;
              cgConnVirtualTmp.push_back(newLid);
              cgConnVirtualGlobal1Tmp.push_back(loc2Glob[oldLid]);
              ++itr;
    
              oldLid = *itr;
              newLid = std::distance(glob2Loc.begin(), glob2Loc.find(loc2Glob[*itr])) + 1;
              old2New[oldLid] = newLid;
              new2Old[newLid] = oldLid;
              cgConnVirtualTmp.push_back(newLid);
              cgConnVirtualGlobal2Tmp.push_back(loc2Glob[oldLid]);
              ++itr;
    
              oldLid = *itr;
              newLid = std::distance(glob2Loc.begin(), glob2Loc.find(loc2Glob[*itr])) + 1;
              old2New[oldLid] = newLid;
              new2Old[newLid] = oldLid;
              cgConnVirtualTmp.push_back(newLid);
              cgConnVirtualGlobal3Tmp.push_back(loc2Glob[oldLid]);
            }
          }
          coords = coordsTmp;
          cgConnReal = cgConnRealTmp;
          cgConnRealGlobal1 = cgConnRealGlobal1Tmp;
          cgConnRealGlobal2 = cgConnRealGlobal2Tmp;
          cgConnRealGlobal3 = cgConnRealGlobal3Tmp;
          if (numVirtualCells)
          {
          cgConnVirtual = cgConnVirtualTmp;
          cgConnVirtualGlobal1 = cgConnVirtualGlobal1Tmp;
          cgConnVirtualGlobal2 = cgConnVirtualGlobal2Tmp;
          cgConnVirtualGlobal3 = cgConnVirtualGlobal3Tmp;
          }
        }

        int numRealPointsInVirtualPatch = 0;
        for (vtkIdType i = it->second->getNumberOfPoints(); i < patchOfPartitionWithVirtualMesh->getNumberOfPoints(); i++)
        {
            patchOfPartitionWithVirtualMesh->getDataSet()->GetPoint(i, pntCoor);
            foundId = pointLocator->FindClosestPoint(pntCoor);
            if ((foundId+1) < partitions[me]->getNumberOfPoints())
            {
              numRealPointsInVirtualPatch++;
            }
        }
        if (!(prefix == "burn" || prefix == "iburn_all"))
        {
          // set num virtual vertices for this partition
          writer->setCoordRind(patchOfPartitionWithVirtualMesh->getNumberOfPoints() // 330 - 326 = 4, total - real = virtual
                               - it->second->getNumberOfPoints() - numRealPointsInVirtualPatch);
          // update number of real vertices in zone:
          writer->setNVrtx(it->second->getNumberOfPoints() + numRealPointsInVirtualPatch);
          writer->setGridXYZ(coords[0], coords[1], coords[2]);
          // Note on why we do the above:
          // The "real" are determined by looking at the real patch. However, it's possible that the virtual faces on the patch have nodes
          // that are "real" on other patches of the same proc. If this is the case, we need to add them to the virtual.
          // Do this by looping through the virtual nodes (last n of the patchOfPartitionWithVirtualMesh, and check if any of them are
          // present in the partitions[me] (local real processor nodes) using the same sort of check as we do below. If this is the case,
          // then we simply accumulate this number. Then subtract it from the rind. The coordinates will still be ordered by global ids,
          // so they are implicitly ordered [real][virtual] according to the global nodes. The element connectivities will be fine, because
          // they reference the order according to global id. Then, after the maps are constructed, the node data will also be assigned
          // the same ids according to the maps. 
        }
        else
        {
          // Remove virtual nodes from coordinates if Rocburn surface
          coords[0].resize(it->second->getNumberOfPoints() + numRealPointsInVirtualPatch);
          coords[1].resize(it->second->getNumberOfPoints() + numRealPointsInVirtualPatch);
          coords[2].resize(it->second->getNumberOfPoints() + numRealPointsInVirtualPatch);
          writer->setGridXYZ(coords[0], coords[1], coords[2]);
          writer->setNVrtx(it->second->getNumberOfPoints() + numRealPointsInVirtualPatch);
        }

        writer->setSection(":t3:real", TRI_3, cgConnReal);

        if (!(prefix == "burn" || prefix == "iburn_all"))
        {
          // get pane data that is stored in element data (patchNo , bcflag, cnstr_type)
          vtkSmartPointer<vtkDataArray> patchNos = 
            patchOfPartitionWithRealMesh->getDataSet()->GetCellData()->GetArray("patchNo");
          double patchNo[1];
          patchNos->GetTuple(0, patchNo);
  
          // Write real and virtual patch information to Rocstar dimension file
          this->dimSurfWriter(me+1, cgConnReal, cgConnVirtual, patchNo[0]);

          vtkSmartPointer<vtkDataArray> bcflags = 
            patchOfPartitionWithRealMesh->getDataSet()->GetCellData()->GetArray("bcflag");
          double bcflag[1];
          bcflags->GetTuple(0, bcflag);
  
          vtkSmartPointer<vtkDataArray> cnstr_types = 
            patchOfPartitionWithRealMesh->getDataSet()->GetCellData()->GetArray("cnstr_type");
          double cnstr_type[1];
          cnstr_types->GetTuple(0, cnstr_type); 

          if (numVirtualCells)
          {
            writer->setNCell(numVirtualCells);
            writer->setSection(":t3:virtual", TRI_3, cgConnVirtual);
            writer->setVirtElmRind(numVirtualCells);
          }
  
          // Before writing Pconn vector, restructure format of vector
          this->restructurePconn(this->surfPconns[me][it->first], me, 1, old2New);
          writer->setPconnVec(this->surfPconns[me][it->first]);
          writer->setPconnLimits(this->pconnProcMin[me], this->pconnProcMax[me]);

          // Note: only tri elements supported
          writer->setGlobalSection("t3g:real#1of3", TRI_3, cgConnRealGlobal1);
          writer->setGlobalNCell(cgConnRealGlobal1.size());
          writer->setGlobalSection("t3g:real#2of3", TRI_3, cgConnRealGlobal2);
          writer->setGlobalNCell(cgConnRealGlobal2.size());
          writer->setGlobalSection("t3g:real#3of3", TRI_3, cgConnRealGlobal3);
          writer->setGlobalNCell(cgConnRealGlobal3.size());
          writer->setGlobalSection("t3g:virtual#1of3", TRI_3, cgConnVirtualGlobal1);
          writer->setGlobalNCell(cgConnVirtualGlobal1.size());
          writer->setGlobalSection("t3g:virtual#2of3", TRI_3, cgConnVirtualGlobal2);
          writer->setGlobalNCell(cgConnVirtualGlobal2.size());
          writer->setGlobalSection("t3g:virtual#3of3", TRI_3, cgConnVirtualGlobal3);
          writer->setGlobalNCell(cgConnVirtualGlobal3.size());

          writer->setGlobalSection("q4g:real#1of4", TETRA_4);
          writer->setGlobalNCell(0);
          writer->setGlobalSection("q4g:real#2of4", TETRA_4);
          writer->setGlobalNCell(0);
          writer->setGlobalSection("q4g:real#3of4", TETRA_4);
          writer->setGlobalNCell(0);
          writer->setGlobalSection("q4g:real#4of4", TETRA_4);
          writer->setGlobalNCell(0);
          writer->setGlobalSection("q4g:virtual#1of4", TETRA_4);
          writer->setGlobalNCell(0);
          writer->setGlobalSection("q4g:virtual#2of4", TETRA_4);
          writer->setGlobalNCell(0);
          writer->setGlobalSection("q4g:virtual#3of4", TETRA_4);
          writer->setGlobalNCell(0);
          writer->setGlobalSection("q4g:virtual#4of4", TETRA_4);
          writer->setGlobalNCell(0);
  
          writer->setPatchNo((int) patchNo[0]);
          writer->setBcflag((int) bcflag[0]);
          writer->setCnstrtype((int) cnstr_type[0]);
        }
        else
        {
          writer->setNCell(1);
          writer->setSection("Empty:t3:virtual", TRI_3, cgConnVirtual);
          writer->setVirtElmRind(numVirtualCells);
        }

        if (!(prefix == "burn" || prefix == "iburn_all"))
        {
          // Subtract surface faces off of volume faces total
          this->nUniqueVolFaces[me] -= patchOfPartitionWithVirtualMesh->getDataSet()->GetNumberOfCells();
        }

        // If the file only contains 1 section, it is empty
        // If we moved this condition furthur up in the file, before all setSection()
        // calls, then we could make it getNZones != 0
        if (!(prefix == "burn" || prefix == "iburn_all"))
        {
          // rocflu surface
          writer->setTypeFlag(1);
          writer->writeZoneToFile();
        }
        else
        {
          // rocburn surface
          writer->setTypeFlag(2);
          writer->writeZoneToFile();
        }

        // write solution data
        if (this->surfWithSol)
        {
          // transfer data to the surf-partition-with-virtual-cells mesh 
          this->burnSurfWithSol
                ->transfer(patchOfPartitionWithRealMesh.get(), "Consistent Interpolation");
          // transfer data to the surf-partition-with-virtual-cells mesh 
          this->surfWithSol
                ->transfer(patchOfPartitionWithVirtualMesh.get(), "Consistent Interpolation");
          int numPointDataArr
            = this->surfWithSol->getDataSet()->GetPointData()->GetNumberOfArrays();
          // begin with point data
          if (numPointDataArr)
          {
            std::string nodeName("NodeData");
            nodeName += this->trimmed_base_t;
            if (!(prefix == "burn" || prefix == "iburn_all"))
            {
              writer->setTypeFlag(0);
              writer->writeSolutionNode(nodeName, Vertex, 0, 1);
              for (int i = 0; i < numPointDataArr; ++i)
              {
                std::vector<double> pointData;
                patchOfPartitionWithVirtualMesh->getPointDataArray(i, pointData);
                std::string dataName(this->surfWithSol->getDataSet()->GetPointData()->GetArrayName(i));
                writer->setTypeFlag(1);
                // Map point data from old coordinates to new coordinates
                this->mapOld2NewPointData(pointData, new2Old);
                writer->writeSolutionField(dataName, nodeName, RealDouble, &pointData[0u]);
              }
            }
            else
            {
              writer->setTypeFlag(0);
              writer->writeSolutionNode(nodeName, Vertex, 1, 1);
            }
          }
          int numCellDataArr;
          if (prefix == "burn" || prefix == "iburn_all")
          {
          numCellDataArr
            = this->burnSurfWithSol->getDataSet()->GetCellData()->GetNumberOfArrays();
          }
          else
          {
          numCellDataArr
            = this->surfWithSol->getDataSet()->GetCellData()->GetNumberOfArrays();            
          }
          if (numCellDataArr)
          {
            if (prefix != "burn")
            {
              std::string nodeName("ElemData");
              nodeName += this->trimmed_base_t;
              if (!(prefix == "burn" || prefix == "iburn_all"))
              {
                writer->setTypeFlag(0);
                writer->writeSolutionNode(nodeName, CellCenter, 0, 1);
              }
              else if (prefix == "iburn_all")
              {
                writer->setTypeFlag(2);
                writer->writeSolutionNode(nodeName, CellCenter, 0, 0);
              }
              // now do cell data
              for (int i = 0; i < numCellDataArr; ++i)
              {
                std::vector<double> cellData;
                std::string dataName;
                if (!(prefix == "burn" || prefix == "iburn_all"))
                {
                  patchOfPartitionWithVirtualMesh->getCellDataArray(i, cellData);
                  dataName = this->surfWithSol->getDataSet()->GetCellData()->GetArrayName(i);
                }
                else if (prefix == "iburn_all")
                {
                  patchOfPartitionWithRealMesh->getCellDataArray(i, cellData);
                  dataName = this->burnSurfWithSol->getDataSet()->GetCellData()->GetArrayName(i);
                }
                if (dataName != "bcflag" && dataName != "patchNo" && dataName != "cnstr_type")
                {
                  // Zero out grid speed for re-meshing
                  if (dataName == "gs")
                  {
                    for (auto itr = cellData.begin(); itr != cellData.end(); ++itr)
                    {
                      *itr = 0;
                    }
                  }

                  // If writing out cell centers, recompute them after re-meshing
                  if (dataName == "centersX" || dataName == "centersY" || dataName == "centersZ")
                  {
                    int ind;
                    if (dataName == "centersX")
                    {
                      ind = 0;
                    }
                    else if (dataName == "centersY")
                    {
                      ind = 1;
                    }
                    else if (dataName == "centersZ")
                    {
                      ind = 2;
                    }
                    std::vector<double> centers;
                    // Get total number of real nodes
                    int nCells = cgConnReal.size()/3;
                    for (int i=0; i<nCells; i++)
                    {
                      centers.push_back((coords[ind][cgConnReal[i]] + coords[ind][cgConnReal[i+1]] + coords[ind][cgConnReal[i+2]])/3);
                    }
                    cellData = centers;
                  }
                  // set num real cells for this patch of this partition
                  if (prefix == "ifluid_ni" || prefix == "ifluid_b")
                  {
                    writer->setTypeFlag(1);
                    if (dataName == "bflag")
                    {
                      writer->writeSolutionField(dataName, nodeName, Integer, &cellData[0u]);
                    }
                    else
                    {
                      writer->writeSolutionField(dataName, nodeName, RealDouble, &cellData[0u]);
                    }
                  }
                  else if (prefix == "burn" || prefix == "iburn_all")
                  {
                    writer->setTypeFlag(2);
                    if (dataName == "bflag")
                    {
                      writer->writeSolutionField(dataName, nodeName, Integer, &cellData[0u]);
                    }
                    else
                    {
                      writer->writeSolutionField(dataName, nodeName, RealDouble, &cellData[0u]);
                    }
                  }
                }
              }
            }
          }
        }
        // reset the section info
        writer->resetSections();
        writer->resetGlobalSections();
      }
    }
    ++it;
  }
  if (patchesFound == 0)
  {
    writer->deleteFile();
  }
}

void RocPartCommGenDriver::writeVolCgns(const std::string& prefix, int proc, int type, int gsExists)
{
  std::stringstream ss;
  ss << prefix << "_" << this->base_t << 
    (proc < 1000 ? (proc < 100 ? (proc < 10 ? "_000" : "_00") : "_0") : "_") << proc << ".cgns";
  std::string cgFname(ss.str());
  std::unique_ptr<cgnsWriter> writer 
    = std::unique_ptr<cgnsWriter>(new cgnsWriter(cgFname, "fluid", 3, 3));
  // initialize Rocstar units, dimensions, and magnitude writing
  writer->setFluidUnitsMap();
  writer->setFluidDimMap();
  writer->setFluidMagMap();
  writer->setiFluidUnitsMap();
  writer->setiFluidDimMap();
  writer->setiFluidMagMap();
  writer->setBurnUnitsMap();
  writer->setBurnDimMap();
  writer->setBurnMagMap();
  // set time stamp
  writer->setTimestamp(this->trimmed_base_t);
  // define elementary information 
  writer->setUnits(Kilogram, Meter, Second, Kelvin, Degree);
  // baseitrname, nTstep, timeVal 
  writer->setBaseItrData("TimeIterValues", 1, std::stod(this->trimmed_base_t));
  // setting zone iter data
  std::string gridpntr("Grid");
  gridpntr += this->trimmed_base_t;
  std::string flowsolpntr("NodeData");
  flowsolpntr += this->trimmed_base_t;
  writer->setZoneItrData("ZoneIterativeData", gridpntr, flowsolpntr);
  ss.str("");
  ss.clear();
  ss << ((proc+1) >= 10 ? "" : "0") << proc+1 << (type < 10 ? "0" : "") << type;
  volZoneMap.push_back(std::stoi(ss.str()));
  // vector to hold all virtual meshes and the partition's mesh for merging
  std::vector<std::shared_ptr<meshBase>> partitionWithAllVirtualCells;
  partitionWithAllVirtualCells.push_back(this->partitions[proc]);
  double myPoint[3];
  for (vtkIdType i = 0; i < partitionWithAllVirtualCells[0]->getDataSet()->GetNumberOfPoints(); i++)
  {
    partitionWithAllVirtualCells[0]->getDataSet()->GetPoint(i, myPoint);
  }

  // iterate over vitual meshes to merge them and get number of virtual coords
  auto virtItr = this->virtualCellsOfPartitions[proc].begin();
  int numVirtualCells = 0;
  while (virtItr != this->virtualCellsOfPartitions[proc].end())
  {
    // push virtual mesh into real
    partitionWithAllVirtualCells.push_back(virtItr->second); 
    numVirtualCells += virtItr->second->getNumberOfCells();
    ++virtItr;
  }

  // merge virtual cells into real mesh for this partition
  std::shared_ptr<meshBase> partitionWithVirtualMesh
    = meshBase::stitchMB(partitionWithAllVirtualCells);

  for (vtkIdType i = 0; i < partitionWithVirtualMesh->getDataSet()->GetNumberOfPoints(); i++)
  {
    partitionWithVirtualMesh->getDataSet()->GetPoint(i, myPoint);
  }

  // define coordinates
  std::vector<std::vector<double>> coords(partitionWithVirtualMesh->getVertCrds());
  // Save vertex coordinates for mapping between volume and surface
  this->procVolMesh.push_back(partitionWithVirtualMesh);
  // set the zone
  writer->setZone(ss.str(), Unstructured);
  // get number of ghost entities in pconn vec
  int ghostDescriptor = this->volPconns[proc].size() - this->notGhostInPconn[proc];
  // register this with writer
  writer->setPconnGhostDescriptor(ghostDescriptor);
  // set num real vertices for this partition
  writer->setNVrtx(partitions[proc]->getNumberOfPoints());
  // set num real cells for this partition
  writer->setNCell(partitions[proc]->getNumberOfCells());
  writer->setGridXYZ(coords[0], coords[1], coords[2]);
  // set num virtual vertices for this partition
  writer->setCoordRind(partitionWithVirtualMesh->getNumberOfPoints()
                       - partitions[proc]->getNumberOfPoints());
  // define real connectivities and 1-base indexing
  std::vector<int> cgConnReal(partitions[proc]->getConnectivities());
  for (auto it = cgConnReal.begin(); it != cgConnReal.end(); ++it)
  {
    *it += 1;
  }
  // define virtual connectivities and 1-base indexing
  std::vector<int> cgConnVirtual(partitionWithVirtualMesh->getConnectivities());
  cgConnVirtual.erase(cgConnVirtual.begin(),cgConnVirtual.begin() + cgConnReal.size());
  for (auto it = cgConnVirtual.begin(); it != cgConnVirtual.end(); ++it)
  {
    *it += 1;
  }
  writer->setSection(":T4:real", TETRA_4, cgConnReal);
  writer->setNCell(numVirtualCells);
  writer->setSection(":T4:virtual", TETRA_4, cgConnVirtual);
  writer->setVirtElmRind(numVirtualCells);

  // Before writing Pconn vector, restructure format of vector
  this->restructurePconn(volPconns[proc], proc, 0, std::map<int,int>());

  writer->setPconnVec(this->volPconns[proc]);
  writer->setPconnLimits(this->pconnProcMin[proc], this->pconnProcMax[proc]);
  writer->setNCell(numVirtualCells);

  if (!gsExists)
  {
    this->dimWriter(proc+1, this->partitions[proc], partitionWithVirtualMesh);
    this->clearBorderVector();
    this->cmpWriter(proc+1, partitionWithVirtualMesh->getDataSet()->GetNumberOfCells());
  }

  // Get number of unique faces in the volume mesh
  if (!gsExists)
  {
    std::vector<std::vector<int>> allFaces;
    std::vector<int> idList;
    for (int iCell = 0; iCell < partitionWithVirtualMesh->getDataSet()->GetNumberOfCells(); iCell++)
    {
      auto myCell = partitionWithVirtualMesh->getDataSet()->GetCell(iCell);
      for (int iFace = 0; iFace < myCell->GetNumberOfFaces(); iFace++)
      {
        auto myFace = myCell->GetFace(iFace);
        idList.clear();
        auto ptIds = myFace->GetPointIds();
        for (int iPt = 0; iPt < myFace->GetNumberOfPoints(); iPt++)
        {
          idList.push_back(ptIds->GetId(iPt));
        }
        std::sort(idList.begin(), idList.end());
        allFaces.push_back(idList);
      }
    }
    std::sort(allFaces.begin(), allFaces.end());
    allFaces.erase(std::unique(allFaces.begin(), allFaces.end()), allFaces.end());
    this->nUniqueVolFaces[proc] = allFaces.size();
  }

  writer->setVolCellFacesNumber(this->nUniqueVolFaces[proc]);

  writer->writeGridToFile();
  writer->setTypeFlag(0);
  writer->writeZoneToFile();

  // write solution data  
  if (this->volWithSol)
  {
    // transfer data to the partition-with-virtual-cells mesh
    this->volWithSol->transfer(partitionWithVirtualMesh.get(),"Consistent Interpolation");
    int numPointDataArr
      = this->volWithSol->getDataSet()->GetPointData()->GetNumberOfArrays();
    if (numPointDataArr)
    {
      std::string nodeName("NodeData");
      nodeName += this->trimmed_base_t;
      writer->setTypeFlag(1);
      writer->writeSolutionNode(nodeName, Vertex, 0, 1); 
      // begin with point data
      for (int i = 0; i < numPointDataArr; ++i)
      {
        std::vector<double> pointData;
        partitionWithVirtualMesh->getPointDataArray(i, pointData);
        std::string dataName(this->volWithSol->getDataSet()->GetPointData()->GetArrayName(i));
        writer->setTypeFlag(0);
        writer->writeSolutionField(dataName, nodeName, RealDouble,&pointData[0u]);
      }
    }
    int numCellDataArr
      = this->volWithSol->getDataSet()->GetCellData()->GetNumberOfArrays();
    if (numCellDataArr)
    {
      std::string nodeName("ElemData"); 
      nodeName += this->trimmed_base_t;
      writer->setTypeFlag(1);
      writer->writeSolutionNode(nodeName, CellCenter, 0, 1); 
      // now do cell data
      for (int i = 0; i < numCellDataArr; ++i)
      {
        std::vector<double> cellData;
        partitionWithVirtualMesh->getCellDataArray(i, cellData);
        std::string dataName(this->volWithSol->getDataSet()->GetCellData()->GetArrayName(i));
        writer->setTypeFlag(0);
        writer->writeSolutionField(dataName, nodeName, RealDouble,&cellData[0u]);
      }
    }
  }
}

void RocPartCommGenDriver::mapOld2NewPointData(std::vector<double> &nodeData, const std::map<int,int> &new2Old)
{
  std::vector<double> nodeDataTmp;
  for (auto itr = new2Old.begin(); itr != new2Old.end(); ++itr)
  {
    nodeDataTmp.push_back(nodeData[itr->second-1]);
  }
  nodeData = nodeDataTmp;
}

void RocPartCommGenDriver::restructurePconn(std::vector<int> &pConnVec, int proc, int volOrSurf, std::map<int,int> old2New)
{
  // Vectors to store each block for each zone
  //<zone, <block, <indices>>
  std::map<int,std::vector<std::vector<int>>> zoneToBlocks;

  // Temporary storage data structures
  std::vector<int> idTmp;

  // Get maximum of PconnVec (maximum of shared nodes)
  int zoneCount = 0;
  int idCount = -1;
  int currZone = -1;
  int max = 0;
  int min = 99999;
  int range = 0;
  int addflag = 0;
  std::vector<int> foundZones;

  for (auto itr = pConnVec.begin(); itr != pConnVec.end(); ++itr)
  {
    if (idCount != -1)
    {
      if (*itr > max && addflag)
      {
        max = *itr;
        pconnProcMax[proc] = max;
      }
      else if (*itr < min && addflag)
      {
        min = *itr;
        pconnProcMin[proc] = min;
      }
      if (!old2New.empty())
      {
        idTmp.push_back(old2New[*itr]);
      }
      else
      {
        idTmp.push_back(*itr);
      }
      idCount++;
    }
    if (zoneCount == 1)
    {
      currZone = *itr;
      auto it = find(foundZones.begin(), foundZones.end(), *itr);
      if (it != foundZones.end())
      {
        addflag = 0;
      }
      else
      {
        addflag = 1;
        foundZones.push_back(*itr);
      }
    }
    if (zoneCount == 2)
    {
      range = *itr;
      idCount = 0;
    }
    if (idCount == range)
    {
      idCount = -1;
      zoneCount =  0;
      zoneToBlocks[currZone].push_back(idTmp);
      idTmp.clear();
    }
    else
    {
      zoneCount++;
    }
  }

  if (volOrSurf == 0)
  {
    std::vector<int> volPconnsTmpProc;
    // Create new Pconn vector
    for (int iBlock = 0; iBlock < 5; iBlock++)
    {  
      volPconnsTmpProc.push_back(zoneToBlocks.size());
      for (auto itr1 = zoneToBlocks.begin(); itr1 != zoneToBlocks.end(); ++itr1)
      {
        volPconnsTmpProc.push_back(itr1->first);
        volPconnsTmpProc.push_back((itr1->second)[iBlock].size());
        for (auto itr2 = (itr1->second)[iBlock].begin(); itr2 != (itr1->second)[iBlock].end(); ++itr2)
        {
          volPconnsTmpProc.push_back(*itr2);
        }
      }
    } 
  pConnVec = volPconnsTmpProc;
  }
  else
  {
    std::vector<int> surfPconnsTmpProc;
    // Create new Pconn vector
    for (int iBlock = 0; iBlock < 1; iBlock++)
    {  
      surfPconnsTmpProc.push_back(zoneToBlocks.size());
      for (auto itr1 = zoneToBlocks.begin(); itr1 != zoneToBlocks.end(); ++itr1)
      {
        surfPconnsTmpProc.push_back(itr1->first);
        surfPconnsTmpProc.push_back((itr1->second)[iBlock].size());
        for (auto itr2 = (itr1->second)[iBlock].begin(); itr2 != (itr1->second)[iBlock].end(); ++itr2)
        {
          surfPconnsTmpProc.push_back(*itr2);
        }
      }
    } 
  pConnVec = surfPconnsTmpProc;
  }
}

void RocPartCommGenDriver::clearPconnVectors()
{
  sharedNodesTmp.clear();
  sentNodesTmp.clear();
  receivedNodesTmp.clear();
  sentCellsTmp.clear();
  receivedCellsTmp.clear();
}

void RocPartCommGenDriver::clearBorderVector()
{
  neighborProcsTmp.clear();
}

// Currently assumes one region per process (i.e. a one-to-one mapping between
// partitions and processes)
void RocPartCommGenDriver::mapWriter()
{
  // Get number of partitions
  int nPart = this->partitions.size();

  // Create mapping file
  std::string fname = this->caseName + ".map";
  ofstream mapFile;
  mapFile.open(fname);

  // Write header
  mapFile << "# ROCFLU region mapping file\n";

  // Write number of regions
  mapFile << "# Number of regions\n";
  mapFile << std::setw(7) << std::to_string(nPart) << std::endl;

  // Write number of processes
  mapFile << "# Number of processes\n";
  mapFile << std::setw(7) << std::to_string(nPart) << std::endl;

  // Write process information
  for (int iPart = 1; iPart < nPart+1; iPart++)
  {
    std::string procStr = std::to_string(iPart);
    std::string procStrPadded = std::string(6 - procStr.length(), '0') + procStr;
    mapFile << "# Process " << procStrPadded << std::endl;

    mapFile << std::setw(7) << std::to_string(1) << std::endl;
    mapFile << std::setw(7) << std::to_string(iPart) << std::endl;
  }

  mapFile << "# End\n";

  mapFile.close();
}

// Cell mapping file
// Currently only supports tetrahedral elements
void RocPartCommGenDriver::cmpWriter(int proc, int nCells)
{
  // Get number of partitions
  int nPart = this->partitions.size();

  // Write cell mapping file
  std::string procStr = std::to_string(proc);
  std::string procStrPadded = std::string(5 - procStr.length(), '0') + procStr;
  std::string fname = this->caseName + ".cmp_" + procStrPadded;
  ofstream cmpFile;
  cmpFile.open(fname);

  // Write header
  cmpFile << "# ROCFLU cell mapping file\n";
  cmpFile << "# Dimensions\n";
  cmpFile << std::setw(8) << std::to_string(nCells)
          << std::setw(8) << std::to_string(0)
          << std::setw(8) << std::to_string(0)
          << std::setw(8) << std::to_string(0) << std::endl;
  cmpFile << "# Tetrahedra\n";
  std::string tmpStr = "";
  std::string cellStr;
  for (int iCell = 1; iCell < nCells+1; iCell++)
  {
    cellStr = std::to_string(iCell);
    tmpStr += std::string(8 - cellStr.length(), ' ') + cellStr;
    if ((iCell) % 10 == 0)
    {
      cmpFile << tmpStr << std::endl;
      tmpStr = "";
    }
  }
  if (tmpStr != "")
  {
    cmpFile << tmpStr << std::endl;
    tmpStr = "";
  }
  cmpFile << "# End" << std::endl;
  cmpFile.close();
}

// Write communication file
// This data is the same as in PaneData/pconn in the fluid CGNS files
void RocPartCommGenDriver::comWriter(int proc)
{

  // Get number of partitions
  int nPart = this->partitions.size();

  // Set number of borders for each proc
  int nBorders[partitions.size()] = {0};
  for (int i = 0; i < partitions.size(); i++)
  {
    nBorders[i] = 1;
  }

  // Write com file for this proc
  std::string procStr = std::to_string(proc);
  std::string procStrPadded = std::string(5 - procStr.length(), '0') + procStr;
  std::string fname = this->caseName + ".com_" + procStrPadded;
  ofstream comFile;
  comFile.open(fname);

  // Write header
  comFile << "# ROCFLU communication lists file\n";

  // Write number of neighboring procs
  comFile << "# Dimensions\n";
  comFile << std::setw(8) << std::to_string(this->neighborProcsTmp.size()) << std::endl;

  // Write number of borders
  comFile << "# Information\n";
  std::string tmpStr;
  std::string borderStr;
  for (auto itr = this->neighborProcsTmp.begin(); itr != neighborProcsTmp.end(); ++itr)
  {
    tmpStr = "";
    procStr = std::to_string((*itr)+1);
    borderStr = std::to_string(nBorders[*itr]);
    tmpStr += std::string(8 - procStr.length(), ' ') + procStr;
    tmpStr += std::string(8 - borderStr.length(), ' ') + borderStr;
    comFile << tmpStr << std::endl;
    nBorders[*itr] += 1;
  }

  // Write cell information
  comFile << "# Cells\n";
  tmpStr = "";
  std::string cellIdStr;
  std::string cellStr1;
  std::string cellStr2;

  //for (auto itr = this->neighborProcsTmp.begin(); itr != neighborProcsTmp.end(); ++itr)
  for (auto itr = this->sentCellsTmp.begin(); itr != sentCellsTmp.end(); ++itr)
  {

    //cellStr1 = std::to_string(this->sentCellsTmp[itr-this->neighborProcsTmp.begin()].size());
    cellStr1 = std::to_string(this->sentCellsTmp[itr-this->sentCellsTmp.begin()].size());
    //cellStr2 = std::to_string(this->receivedCellsTmp[itr-this->neighborProcsTmp.begin()].size());
    cellStr2 = std::to_string(this->receivedCellsTmp[itr-this->sentCellsTmp.begin()].size());
    tmpStr += std::string(8 - cellStr1.length(), ' ') + cellStr1;
    tmpStr += std::string(8 - cellStr2.length(), ' ') + cellStr2;
    comFile << tmpStr << std::endl;
    tmpStr = "";
    for (auto itr2 = this->sentCellsTmp[itr-this->sentCellsTmp.begin()].begin();
        itr2 < this->sentCellsTmp[itr-this->sentCellsTmp.begin()].end();
        ++itr2)
    {
      cellIdStr = std::to_string(*itr2);
      tmpStr += std::string(8 - cellIdStr.length(), ' ') + cellIdStr;
      if (((itr2-this->sentCellsTmp[itr-this->sentCellsTmp.begin()].begin())+1) % 10 == 0)
      {
        comFile << tmpStr << std::endl;
        tmpStr = "";
      }
    }
    if (tmpStr != "")
    {
      comFile << tmpStr << std::endl;
      tmpStr = "";
    }
    for (auto itr2 = this->receivedCellsTmp[itr-this->sentCellsTmp.begin()].begin();
        itr2 < this->receivedCellsTmp[itr-this->sentCellsTmp.begin()].end();
        ++itr2)
    {
      cellIdStr = std::to_string(*itr2);
      tmpStr += std::string(8 - cellIdStr.length(), ' ') + cellIdStr;
      if (((itr2-this->receivedCellsTmp[itr-this->sentCellsTmp.begin()].begin())+1) % 10 == 0)
      {        
        comFile << tmpStr << std::endl;
        tmpStr = "";
      }
    }
    if (tmpStr != "")
    {
      comFile << tmpStr << std::endl;
      tmpStr = "";
    }
  }

  // Write node information
  comFile << "# Vertices\n";
  tmpStr = "";
  std::string vertIdStr;
  std::string vertStr1;
  std::string vertStr2;
  std::string vertStr3;

  for (auto itr = this->neighborProcsTmp.begin(); itr != neighborProcsTmp.end(); ++itr)
  {

    vertStr1 = std::to_string(this->sentNodesTmp[itr-this->neighborProcsTmp.begin()].size());
    vertStr2 = std::to_string(this->receivedNodesTmp[itr-this->neighborProcsTmp.begin()].size());
    vertStr3 = std::to_string(this->sharedNodesTmp[itr-this->neighborProcsTmp.begin()].size());
    tmpStr += std::string(8 - vertStr1.length(), ' ') + vertStr1;
    tmpStr += std::string(8 - vertStr2.length(), ' ') + vertStr2;
    tmpStr += std::string(8 - vertStr3.length(), ' ') + vertStr3;
    comFile << tmpStr << std::endl;
    tmpStr = "";
    for (auto itr2 = this->sentNodesTmp[itr-this->neighborProcsTmp.begin()].begin();
        itr2 < this->sentNodesTmp[itr-this->neighborProcsTmp.begin()].end();
        ++itr2)
    {
      vertIdStr = std::to_string(*itr2);
      tmpStr += std::string(8 - vertIdStr.length(), ' ') + vertIdStr;
      if (((itr2-this->sentNodesTmp[itr-this->neighborProcsTmp.begin()].begin())+1) % 10 == 0)
      {
       comFile << tmpStr << std::endl;
       tmpStr = "";
      }
    }
    if (tmpStr != "")
    {
      comFile << tmpStr << std::endl;
      tmpStr = "";
    }
    for (auto itr2 = this->receivedNodesTmp[itr-this->neighborProcsTmp.begin()].begin();
        itr2 < this->receivedNodesTmp[itr-this->neighborProcsTmp.begin()].end();
        ++itr2)
    {
      vertIdStr = std::to_string(*itr2);
      tmpStr += std::string(8 - vertIdStr.length(), ' ') + vertIdStr;
      if (((itr2-this->receivedNodesTmp[itr-this->neighborProcsTmp.begin()].begin())+1) % 10 == 0)
      {
       comFile << tmpStr << std::endl;
       tmpStr = "";
      }
    }
    if (tmpStr != "")
    {
      comFile << tmpStr << std::endl;
      tmpStr = "";
    }
    for (auto itr2 = this->sharedNodesTmp[itr-this->neighborProcsTmp.begin()].begin();
        itr2 < this->sharedNodesTmp[itr-this->neighborProcsTmp.begin()].end();
        ++itr2)
    {
      vertIdStr = std::to_string(*itr2);
      tmpStr += std::string(8 - vertIdStr.length(), ' ') + vertIdStr;
      if (((itr2-this->sharedNodesTmp[itr-this->neighborProcsTmp.begin()].begin())+1) % 10 == 0)
      {
       comFile << tmpStr << std::endl;
       tmpStr = "";
      }
    }
    if (tmpStr != "")
    {
      comFile << tmpStr << std::endl;
      tmpStr = "";
    }
  }

  comFile << "# End" << std::endl;
  comFile.close();
}

// Rocflu dimensions file
void RocPartCommGenDriver::dimWriter(int proc, std::shared_ptr<meshBase> realMB, std::shared_ptr<meshBase> totalMB)
{

  // Set number of borders for each proc
  int nBorders[partitions.size()] = {0};
  for (int i = 0; i < partitions.size(); i++)
  {
    nBorders[i] = 1;
  }

  // Write com file for this proc
  std::string procStr = std::to_string(proc);
  std::string procStrPadded = std::string(5 - procStr.length(), '0') + procStr;
  std::string fname = this->caseName + ".dim_" + procStrPadded + "_0.00000E+00";
  ofstream dimFile;
  dimFile.open(fname);

  // Write vertices
  dimFile << "# ROCFLU dimensions file\n";

  // Write cells
  dimFile << "# Vertices\n";
  dimFile << std::setw(8) << std::to_string(realMB->getDataSet()->GetNumberOfPoints());
  dimFile << std::setw(8) << std::to_string(totalMB->getDataSet()->GetNumberOfPoints());
  dimFile << std::setw(8) << std::to_string(totalMB->getDataSet()->GetNumberOfPoints()*4) << std::endl;

  // Write header
  dimFile << "# Cells\n";
  dimFile << std::setw(8) << std::to_string(realMB->getDataSet()->GetNumberOfCells());
  dimFile << std::setw(8) << std::to_string(totalMB->getDataSet()->GetNumberOfCells());
  dimFile << std::setw(8) << std::to_string(totalMB->getDataSet()->GetNumberOfCells()*4) << std::endl;

  // Write tetrahedra (currently supports only tetrahedra)
  dimFile << "# Tetrahedra\n";
  dimFile << std::setw(8) << std::to_string(realMB->getDataSet()->GetNumberOfCells());
  dimFile << std::setw(8) << std::to_string(totalMB->getDataSet()->GetNumberOfCells());
  dimFile << std::setw(8) << std::to_string(totalMB->getDataSet()->GetNumberOfCells()*4) << std::endl;

  // Write hexahedra (not supported)
  dimFile << "# Hexahedra\n";
  dimFile << std::setw(8) << "0";
  dimFile << std::setw(8) << "0";
  dimFile << std::setw(8) << "0" << std::endl;

  // Write prisms (not supported)
  dimFile << "# Prisms\n";
  dimFile << std::setw(8) << "0";
  dimFile << std::setw(8) << "0";
  dimFile << std::setw(8) << "0" << std::endl;

  // Write pyramids (not supported)
  dimFile << "# Pyramids\n";
  dimFile << std::setw(8) << "0";
  dimFile << std::setw(8) << "0";
  dimFile << std::setw(8) << "0" << std::endl;

  // Write patches
  dimFile << "# Patches (v2)\n";

  // Write borders
  dimFile << "# Borders\n";
  if (proc == 0)
  {
    dimFile << std::setw(8) << "0" << std::endl;
  }
  else
  {
    // Write number of neighboring procs
    dimFile << std::setw(8) << std::to_string(this->neighborProcsTmp.size()) << std::endl;
    
    // Write number of borders
    std::string tmpStr;
    std::string borderStr;
    std::string cellStr1;
    std::string cellStr2;
    std::string vertStr1;
    std::string vertStr2;
    std::string vertStr3;
    for (auto itr = this->neighborProcsTmp.begin(); itr != neighborProcsTmp.end(); ++itr)
    {
      tmpStr = "";

      // proc/border information
      procStr = std::to_string((*itr)+1);
      borderStr = std::to_string(nBorders[*itr]);
      tmpStr += std::string(8 - procStr.length(), ' ') + procStr;
      tmpStr += std::string(8 - borderStr.length(), ' ') + borderStr;
      nBorders[*itr] += 1;

      // sent, received cells
      cellStr1 = std::to_string(this->sentCellsTmp[itr-this->neighborProcsTmp.begin()].size());
      cellStr2 = std::to_string(this->receivedCellsTmp[itr-this->neighborProcsTmp.begin()].size());
      tmpStr += std::string(8 - cellStr1.length(), ' ') + cellStr1;
      tmpStr += std::string(8 - cellStr2.length(), ' ') + cellStr2;

      // sent, received, shared vertices
      vertStr1 = std::to_string(this->sentNodesTmp[itr-this->neighborProcsTmp.begin()].size());
      vertStr2 = std::to_string(this->receivedNodesTmp[itr-this->neighborProcsTmp.begin()].size());
      vertStr3 = std::to_string(this->sharedNodesTmp[itr-this->neighborProcsTmp.begin()].size());
      tmpStr += std::string(8 - vertStr1.length(), ' ') + vertStr1;
      tmpStr += std::string(8 - vertStr2.length(), ' ') + vertStr2;
      tmpStr += std::string(8 - vertStr3.length(), ' ') + vertStr3;

      dimFile << tmpStr << std::endl;
    }
  }

  // Write borders
  dimFile << "# End\n";
  dimFile.close();
}

// Write surface information for Rocstar dimension file
void RocPartCommGenDriver::dimSurfWriter(int proc, std::vector<int> cgConnReal, std::vector<int> cgConnVirtual, int patchNo)
{
  // Write com file for this proc
  std::string procStr = std::to_string(proc);
  std::string procStrPadded = std::string(5 - procStr.length(), '0') + procStr;
  std::string fname = this->caseName + ".dim_" + procStrPadded + "_0.00000E+00";

  std::string line;
  std::vector< std::string > myLines;
  std::ifstream myfile(fname);

  // Get current dim file
  while (std::getline(myfile, line))
  {
    myLines.push_back(line);
  }

  // Accumulate all patches
  std::set<int> patches;
  for (auto itr = patchesOfSurfacePartitions.begin(); itr != patchesOfSurfacePartitions.end(); ++itr)
  {
    for (auto itr2 = (*itr).begin(); itr2 != (*itr).end(); ++itr2)
    {
      patches.insert(itr2->first);
    }
  }

  // Insert Patch lines in existing file
  std::vector < std::string > newLines;
  std::string patchStr1;
  std::string patchStr2;
  std::string patchStr3;
  std::string patchStr4;
  std::string tmpStr;
  for (auto itr = myLines.begin(); itr != myLines.end(); ++itr)
  {
    //if ((*itr).find("Patches (v2)") != std::string::npos)
    if ((*itr).find("Borders") != std::string::npos)
    {
      if ((*(itr-1)).find("Patches (v2)") != std::string::npos)
      {
        tmpStr = "";
        patchStr1 = std::to_string(this->patchesOfSurfacePartitions[proc-1].size());
        patchStr2 = std::to_string(*patches.rbegin());
        tmpStr += std::string(8 - patchStr1.length(), ' ') + patchStr1;
        tmpStr += std::string(8 - patchStr2.length(), ' ') + patchStr2;
        newLines.push_back(tmpStr);
      }
      for (auto itr2 = this->patchesOfSurfacePartitions[proc-1].begin(); itr2 != this->patchesOfSurfacePartitions[proc-1].end(); ++itr2)
      {
        if (itr2->first == patchNo)
        {
          this->totalTrisPerPatch[patchNo] += cgConnReal.size()+cgConnVirtual.size();
          tmpStr = "";
          patchStr1 = std::to_string(itr2->first);
          tmpStr += std::string(8 - patchStr1.length(), ' ') + patchStr1;
          patchStr2 = std::to_string(cgConnReal.size()/3);
          tmpStr += std::string(8 - patchStr2.length(), ' ') + patchStr2;
          patchStr3 = std::to_string((cgConnReal.size()+cgConnVirtual.size())/3);
          tmpStr += std::string(8 - patchStr3.length(), ' ') + patchStr3;
          patchStr4 = std::to_string(4*(cgConnReal.size()+cgConnVirtual.size())/3);
          tmpStr += std::string(8 - patchStr4.length(), ' ') + patchStr4;
          tmpStr += std::string(8 - 1, ' ') + '0';
          tmpStr += std::string(8 - 1, ' ') + '0';
          tmpStr += std::string(8 - 1, ' ') + '0';
          tmpStr += std::string(8 - 1, ' ') + '0';
          newLines.push_back(tmpStr);
        }
      }
    }
    newLines.push_back(*itr);
  }

  // Write new lines to file
  ofstream dimFile;
  dimFile.open(fname);
  for (auto itr = newLines.begin(); itr != newLines.end(); ++itr)
  {
    dimFile << *itr << std::endl;
  }

  dimFile.close();
}


// Write surface information for Rocstar dimension file
void RocPartCommGenDriver::dimSurfWriter(int proc)
{

  // Write com file for this proc
  std::string procStr = std::to_string(proc);
  std::string procStrPadded = std::string(5 - procStr.length(), '0') + procStr;
  std::string fname = this->caseName + ".dim_" + procStrPadded + "_0.00000E+00";

  std::string line;
  std::vector< std::string > myLines;
  std::ifstream myfile(fname);

  // Get current dim file
  while (std::getline(myfile, line))
  {
    myLines.push_back(line);
  }

  // Accumulate all patches
  std::set<int> patches;
  for (auto itr = patchesOfSurfacePartitions.begin(); itr != patchesOfSurfacePartitions.end(); ++itr)
  {
    for (auto itr2 = (*itr).begin(); itr2 != (*itr).end(); ++itr2)
    {
      patches.insert(itr2->first);
    }
  }

  // Insert Patch lines in existing file
  std::vector < std::string > newLines;
  std::string patchStr1;
  std::string patchStr2;
  std::string patchStr3;
  std::string patchStr4;
  std::string tmpStr;
  for (auto itr = myLines.begin(); itr != myLines.end(); ++itr)
  {
    if ((*itr).find("Borders") != std::string::npos)
    {
      if ((*(itr-1)).find("Patches (v2)") != std::string::npos)
      {
        tmpStr = "";
        patchStr1 = std::to_string(*patches.rbegin());
        tmpStr += std::string(8 - patchStr1.length(), ' ') + patchStr1;
        tmpStr += std::string(8 - patchStr1.length(), ' ') + patchStr1;
        newLines.push_back(tmpStr);
      }
      for (auto itr2 = patches.begin(); itr2 != patches.end(); ++itr2)
      {
        tmpStr = "";
        patchStr1 = std::to_string(*itr2);
        tmpStr += std::string(8 - patchStr1.length(), ' ') + patchStr1;
        patchStr2 = std::to_string((this->totalTrisPerPatch[*itr2])/3);
        tmpStr += std::string(8 - patchStr2.length(), ' ') + patchStr2;
        patchStr3 = std::to_string((this->totalTrisPerPatch[*itr2])/3);
        tmpStr += std::string(8 - patchStr3.length(), ' ') + patchStr3;
        patchStr4 = std::to_string(4*(this->totalTrisPerPatch[*itr2])/3);
        tmpStr += std::string(8 - patchStr4.length(), ' ') + patchStr4;
        tmpStr += std::string(8 - 1, ' ') + '0';
        tmpStr += std::string(8 - 1, ' ') + '0';
        tmpStr += std::string(8 - 1, ' ') + '0';
        tmpStr += std::string(8 - 1, ' ') + '0';
        newLines.push_back(tmpStr);
      }
    }
    newLines.push_back(*itr);
  }

  // Write new lines to file
  ofstream dimFile;
  dimFile.open(fname);
  for (auto itr = newLines.begin(); itr != newLines.end(); ++itr)
  {
    dimFile << *itr << std::endl;
  }

  dimFile.close();

}

void RocPartCommGenDriver::txtWriter()
{
  // Write txt files for fluid and ifluid
  std::string fname;
  fname = "fluid_in_" + this->base_t + ".txt";
  ofstream fluid_in;
  fluid_in.open(fname);

  fname = "ifluid_in_" + this->base_t + ".txt";
  ofstream ifluid_in;
  ifluid_in.open(fname);

  fname = "burn_in_" + this->base_t + ".txt";
  ofstream burn_in;
  burn_in.open(fname);

  fname = "iburn_in_" + this->base_t + ".txt";
  ofstream iburn_in;
  iburn_in.open(fname);


  // Get number of partitions
  int nPart = this->partitions.size();

  std::string tmpStr;
  std::string paneStr;
  for (int iPart = 0; iPart < nPart; iPart++)
  {
    // Write to fluid file
    fluid_in << "@Proc: " << std::to_string(iPart) << std::endl;
    fluid_in << "@Files: fluid_" + this->base_t + "_%04p.cgns" << std::endl;
    tmpStr = "@Panes:";
    tmpStr += " " + std::to_string(volZoneMap[iPart])+"\n";
    fluid_in << tmpStr << std::endl;

    // Write to ifluid_file
    ifluid_in << "@Proc: " << std::to_string(iPart) << std::endl;;
    ifluid_in << "@Files: ifluid*_" + this->base_t + "_%04p.cgns" << std::endl;
    tmpStr = "@Panes:";
    for (auto itr1 = surfZoneMap[iPart].begin(); itr1 != surfZoneMap[iPart].end(); ++itr1)
    {
      tmpStr += " " + std::to_string(iPart+1);
      if ((itr1->second) < 10)
      {
        tmpStr += "0";
      }
      tmpStr += std::to_string(itr1->second);
    }
    tmpStr += "\n";
    ifluid_in << tmpStr << std::endl;


    // Write to iburn_file
    iburn_in << "@Proc: " << std::to_string(iPart) << std::endl;;
    tmpStr = "@Panes:";
    int burnCnt = 0;
    for (auto itr1 = surfZoneMap[iPart].begin(); itr1 != surfZoneMap[iPart].end(); ++itr1)
    {
      if (itr1->first == 1)
      {
        burnCnt += 1;
        iburn_in << "@Files: iburn*_" + this->base_t + "_%04p.cgns" << std::endl;
        tmpStr += " " + std::to_string(iPart+1);
      if ((itr1->second) < 10)
        {
          tmpStr += "0";
        }
        tmpStr += std::to_string(itr1->second);
      }
    }
    if (burnCnt == 0)
    {
      iburn_in << "@Files:" << std::endl;
    }
    tmpStr += "\n";
    iburn_in << tmpStr << std::endl;

    // Write to iburn_file
    burn_in << "@Proc: " << std::to_string(iPart) << std::endl;;
    tmpStr = "@Panes:";
    burnCnt = 0;
    for (auto itr1 = surfZoneMap[iPart].begin(); itr1 != surfZoneMap[iPart].end(); ++itr1)
    {
      if (itr1->first == 1)
      {
        burnCnt += 1;
        burn_in << "@Files: burn_" + this->base_t + "_%04p.cgns" << std::endl;
        tmpStr += " " + std::to_string(iPart+1);
        if ((itr1->second) < 10)
        {
          tmpStr += "0";
        }
        tmpStr += std::to_string(itr1->second);
      }
    }
    if (burnCnt == 0)
    {
      burn_in << "@Files:" << std::endl;
    }
    tmpStr += "\n";
    burn_in << tmpStr << std::endl;
  }
  fluid_in.close();
  ifluid_in.close();
}

void RocPartCommGenDriver::swapTriOrdering(std::vector<int> &connVec)
{
  // Change ordering from counter-clockwise to clockwise convention
  std::vector<int> connVecTmp;
  int ind = 1;
  int tmpInd = -1;
  for (auto itr = connVec.begin(); itr != connVec.end(); ++itr)
  {
    if (ind == 1)
    {
      connVecTmp.push_back(*itr);
      ind++;
    }
    else if (ind == 2)
    {
      tmpInd = *itr;
      ind++;
    }
    else if (ind == 3)
    {
      connVecTmp.push_back(*itr);
      connVecTmp.push_back(tmpInd);
      ind = 1;
    }
  }
  connVec = connVecTmp;
}
