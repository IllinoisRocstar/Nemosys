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

#include <sstream>
#include <cstddef>
#include <AuxiliaryFunctions.H>
#include <cgnsWriter.H>

/*
  TODO: handling the iburn and burn files
        - 



*/

RocPartCommGenDriver::RocPartCommGenDriver(std::shared_ptr<meshBase> _mesh, 
                                           std::shared_ptr<meshBase> _remeshedSurf, 
                                           std::shared_ptr<meshBase> _volWithSol,
                                           std::shared_ptr<meshBase> _surfWithSol,
                                           int numPartitions, const std::string& _base_t,
                                           bool _writeIntermediateFiles,
                                           double _searchTolerance)
{
  std::cout << "RocPartCommGenDriver created\n";
  this->mesh = _mesh;
  // load stitched surf mesh with patch info
  this->remeshedSurf = _remeshedSurf;
  this->volWithSol = _volWithSol;
  this->surfWithSol = _surfWithSol;
  this->base_t = _base_t; 
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
  remeshedSurf->setContBool(false);
  this->partitions = meshBase::partition(this->mesh.get(), numPartitions);
  this->AddGlobalCellIds(this->remeshedSurf);
  std::vector<std::string> patchNoAndGlobalCellIds = {"patchNo","GlobalCellIds"};
  // initialize storage for surf partitions 
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  this->surfacePartitions.resize(numPartitions);
  // initialize storage for virtual cells
  this->virtualCellsOfPartitions.resize(numPartitions);
  this->virtualCellsOfSurfPartitions.resize(numPartitions);
  // get all required local-global node maps and node/cell ids
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
                           "Consistent Interpolation", patchNoAndGlobalCellIds, 1);
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
    this->writeVolCgns("fluid", i,1);
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
  for (int i = 0; i < numPartitions; ++i)
  {
    this->writeSharedToPconn(i);
    this->writeSurfCgns("ifluid_ni", i);
  }
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
      std::stringstream ss;
      ss << sharedItr->first + 1 << type;
      // zone number
      this->volPconns[proc].push_back(std::stoi(ss.str()));
      // number of ids
      this->volPconns[proc].push_back(sharedItr->second.size());
      // indices 
      for (int i = 0; i < sharedItr->second.size(); ++i)
      {
        this->volPconns[proc].push_back(sharedItr->second[i] + 1);
      }
    }
    ++sharedItr;
  }
  return volPconns[proc].size();
}

void RocPartCommGenDriver::writeSharedToPconn(int proc)
{
  auto it = this->sharedPatchNodes[proc].begin();
  while (it != sharedPatchNodes[proc].end())
  {
    auto it1 = it->second.begin();
    while (it1 != it->second.end())
    {
      auto it2 = it1->second.begin();
      while (it2 != it1->second.end())
      {
        // num blocks for this proc and current patch
        this->surfPconns[proc][it->first].push_back(1);
        // zone number
        std::stringstream ss;
        ss << it1->first + 1 << "0" << it2->first;
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
  	  while (tmpIt != sentItr->second.end())
  	  {
        this->volPconns[proc].push_back(*tmpIt + 1);
  	  	++tmpIt;
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
    for (int id = numPrevReceived; id < numReceived; ++id)
    {
      this->volPconns[proc].push_back(offset+id+1);
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
  meshBase* meMesh = (vol ? partitions[me].get() : surfacePartitions[me].get());

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
    // for each cell using shared node (these will be cells on the boundary)
    for (int k = 0; k < cellIdsList->GetNumberOfIds(); ++k)
    {
      // get the local cell idx
      int localCellId = cellIdsList->GetId(k);
      // get the cell in me for point extraction
      meMesh->getDataSet()->GetCell(localCellId, genCell);
      int numSharedInCell = 0;
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
            // add idx to map of nodes me sends to you
            this->sentNodes[me][you].insert(pntId);
            // add idx to map of nodes you recevies from proc me 
            this->receivedNodes[you][me].insert(globPntId);
          }
          else
          {
            // add idx to map of nodes me sends to you
            this->sentSurfNodes[me][you].insert(pntId);
            // add idx to map of nodes you recevies from proc me 
            this->receivedSurfNodes[you][me].insert(globPntId);
          }
        }
        else 
        {
          numSharedInCell += 1;
        }
        // if there are more than two shared node with other proc, we send the cell
        if (numSharedInCell >= 2)
        {
          if (vol)
          {
            // add idx to sentCells to you map
            this->sentCells[me][you].insert(localCellId);
          }
          else
          {
            this->sentSurfCells[me][you].insert(localCellId);
          }
          // get the global cell index
          int globId = partToGlobCellMap[me][localCellId];
          if (vol)
          {
            // add idx to map of cells proc you recieves from proc me
            this->receivedCells[you][me].insert(globId);
          }
          else
          {
            this->receivedSurfCells[you][me].insert(globId);
          }
        }
      }
    }
  }
}

void RocPartCommGenDriver::getSharedPatchInformation()
{
  // get glob to part node maps for intersection calc
  std::map<int, std::vector<std::map<int,int>>> patchProcGlobToPartNodeMaps;
  std::map<int, std::vector<std::map<int,int>>> patchProcPartToGlobNodeMaps;
  std::map<int, std::vector<std::vector<int>>> patchProcGlobNodeIdsMaps;

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
      patchProcGlobToPartNodeMaps[it->first].push_back(patchGlobToPartNodeMap);
      patchProcPartToGlobNodeMaps[it->first].push_back(patchPartToGlobNodeMap);
      patchProcGlobNodeIdsMaps[it->first].push_back(getSortedKeys(patchGlobToPartNodeMap)); 
      ++it;
    }
  }
  for (int me = 0; me < this->surfacePartitions.size(); ++me)
  {
    auto it = this->patchesOfSurfacePartitions[me].begin();
    while (it != this->patchesOfSurfacePartitions[me].end())
    {   
      if (me < patchProcGlobNodeIdsMaps[it->first].size())
      {  
        for (int you = 0; you < this->surfacePartitions.size(); ++you)
        {
          auto it1 = this->patchesOfSurfacePartitions[you].begin();

          while (it1 != this->patchesOfSurfacePartitions[you].end())
          {
            if (you < patchProcGlobNodeIdsMaps[it1->first].size())   
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
                    // get local id in patch in proc me
                    int localPntId = patchProcGlobToPartNodeMaps[it->first][me][tmpVec[i]];
                    this->sharedPatchNodes[me][it->first][you][it1->first][i] = localPntId;
                    // get local id in patch in proc you
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
        std::cout << "virtual patch num " << it1->first << " of proc " << i << " from proc " << it->first
                  << "\n\n" << std::endl;
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
  // define elementary information
  writer->setUnits(Kilogram, Meter, Second, Kelvin, Degree);
  // baseitrname, nTstep, timeVal
  writer->setBaseItrData("TimeIterValues", 1, std::stod(this->trimmed_base_t));
  // write elementary grid information to file
  writer->writeGridToFile();
  std::string grdpntr("Grid");
  grdpntr += this->trimmed_base_t;
  std::string flowsolpntr("Grid");
  flowsolpntr += this->trimmed_base_t;
  ss.str("");
  ss.clear();
  //TODO : // determine whether this surface partition has a burning patch
  int zone = 2;
  if (this->patchesOfSurfacePartitions[me].find(1) 
        != this->patchesOfSurfacePartitions[me].end())
  {
    zone += 1; // zone starts at 3 in case there is a burning patch on partition 
  }
  // begin loop over patches of this partition
  auto it = this->patchesOfSurfacePartitions[me].begin();
  while (it != this->patchesOfSurfacePartitions[me].end())
  {
    if (it->first != 1) // ni files have zones > 1 and 1 is a burning patch number
    {
      int numVirtualCells = 0;;
      // setting zone iter data
      writer->setZoneItrData("ZoneIterativeData", grdpntr, flowsolpntr);
      ss.str("");
      ss.clear();
      ss << 0 << me+1 << (it->first < 10 ? "0" : "") <<  zone;
      zone += 1;    
      // vector to hold real patch mesh and all virtuals from same patch on other procs
      std::vector<std::shared_ptr<meshBase>> patchOfPartitionWithAllVirtualCells;
      // first push back real cells
      patchOfPartitionWithAllVirtualCells.push_back(it->second);
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
      // patchOfPartitionWithAllVirtualCells
      // first meshBase is real: "extractedPatch*.vtu"
      // next n meshBases are virtual: "extractedVirtualPatch*.vtu"
      if (numVirtualCells)
      { 
        // merge virtual cells into real mesh for this partition
        std::shared_ptr<meshBase> patchOfPartitionWithVirtualMesh
          = meshBase::stitchMB(patchOfPartitionWithAllVirtualCells);
        // set the zone
        writer->setZone(ss.str(), Unstructured);
        // all pconn stuff is shared nodes, so no ghost entities
        writer->setPconnGhostDescriptor(0);
        // set num real vertices
        writer->setNVrtx(it->second->getNumberOfPoints());
        // set num real cells for this patch of this partition
        writer->setNCell(it->second->getNumberOfCells());
        // define and set coordinates
        std::vector<std::vector<double>> 
          coords(patchOfPartitionWithVirtualMesh->getVertCrds());
        writer->setGridXYZ(coords[0], coords[1], coords[2]);
        // set num virtual vertices for this partition
        writer->setCoordRind(patchOfPartitionWithVirtualMesh->getNumberOfPoints()
                             - it->second->getNumberOfPoints());
        // define real connectivities and 1-base indexing
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
        writer->setSection(":T3:real", TRI_3, cgConnReal);
        writer->setNCell(numVirtualCells);
        writer->setSection(":T3:virtual", TRI_3, cgConnVirtual);
        writer->setVirtElmRind(numVirtualCells);
        writer->setPconnVec(this->surfPconns[me][it->first]);
        writer->writeZoneToFile();
        // write solution data
        if (this->surfWithSol)
        {   
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
            writer->writeSolutionNode(nodeName, Vertex);
            for (int i = 0; i < numPointDataArr; ++i)
            {
              std::vector<double> pointData;
              patchOfPartitionWithVirtualMesh->getPointDataArray(i, pointData);
              std::string dataName(this->surfWithSol->getDataSet()->GetPointData()->GetArrayName(i));
              writer->writeSolutionField(dataName, nodeName, RealDouble, &pointData[0u]);
            } 
          }
          int numCellDataArr
            = this->surfWithSol->getDataSet()->GetCellData()->GetNumberOfArrays();
          if (numCellDataArr)
          {
            std::string nodeName("ElemData"); 
            nodeName += this->trimmed_base_t;
            writer->writeSolutionNode(nodeName, CellCenter); 
            // now do cell data
            for (int i = 0; i < numCellDataArr; ++i)
            {
              std::vector<double> cellData;
              patchOfPartitionWithVirtualMesh->getCellDataArray(i, cellData);
              std::string dataName(this->surfWithSol->getDataSet()->GetCellData()->GetArrayName(i));
              // TODO: these fields are stored as cell data in the stitchedSurf.vtu, i.e. surfWithSol
              //       we don't want to write it back to cgns as solution data. instead, grab one of the 
              //       values from each using something like:
    // patchOfPartitionWithVirtualMesh->getDataSet()->GetPointData()->GetArray("patchNo")->GetTuple(0);
              // add a line like above BEFORE writeZoneToFile
              // use set methods and add members to cgnsWriter and add the data to PaneData in the 
              // writeZoneToFile function of cgnsWriter (see PaneData example there)
              if (dataName != "bcflag" && dataName != "patchNo" && dataName != "cnstr_type")
              {
                writer->writeSolutionField(dataName, nodeName, RealDouble,&cellData[0u]);
              }
            }
          }
        }
        // reset the section info
        writer->resetSections();
      }
    }
    ++it;
  }
}

void RocPartCommGenDriver::writeVolCgns(const std::string& prefix, int proc, int type)
{
  std::stringstream ss;
  ss << prefix << "_" << this->base_t << 
    (proc < 1000 ? (proc < 100 ? (proc < 10 ? "_000" : "_00") : "_0") : "_") << proc << ".cgns";

  std::string cgFname(ss.str());
  std::unique_ptr<cgnsWriter> writer 
    = std::unique_ptr<cgnsWriter>(new cgnsWriter(cgFname, "fluid", 3, 3));
  // define elementary information 
  writer->setUnits(Kilogram, Meter, Second, Kelvin, Degree);
  // baseitrname, nTstep, timeVal 
  writer->setBaseItrData("TimeIterValues", 1, std::stod(this->trimmed_base_t));
  // setting zone iter data
  std::string gridpntr("Grid");
  gridpntr += this->trimmed_base_t;
  std::string flowsolpntr("Grid");
  flowsolpntr += this->trimmed_base_t;
  writer->setZoneItrData("ZoneIterativeData", gridpntr, flowsolpntr);
  ss.str("");
  ss.clear();
  ss << 0 << proc+1 << (type < 10 ? "0" : "") << type;
  // vector to hold all virtual meshes and the partition's mesh for merging
  std::vector<std::shared_ptr<meshBase>> partitionWithAllVirtualCells;
  partitionWithAllVirtualCells.push_back(this->partitions[proc]);
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
  // define coordinates
  std::vector<std::vector<double>> coords(partitionWithVirtualMesh->getVertCrds());
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
  writer->setPconnVec(this->volPconns[proc]);
  writer->writeGridToFile();
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
      writer->writeSolutionNode(nodeName, Vertex); 
      // begin with point data
      for (int i = 0; i < numPointDataArr; ++i)
      {
        std::vector<double> pointData;
        partitionWithVirtualMesh->getPointDataArray(i, pointData);
        std::string dataName(this->volWithSol->getDataSet()->GetPointData()->GetArrayName(i));
        writer->writeSolutionField(dataName, nodeName, RealDouble,&pointData[0u]);
      }
    }
    int numCellDataArr
      = this->volWithSol->getDataSet()->GetCellData()->GetNumberOfArrays();
    if (numCellDataArr)
    {
      std::string nodeName("ElemData"); 
      nodeName += this->trimmed_base_t;
      writer->writeSolutionNode(nodeName, CellCenter); 
      // now do cell data
      for (int i = 0; i < numCellDataArr; ++i)
      {
        std::vector<double> cellData;
        partitionWithVirtualMesh->getCellDataArray(i, cellData);
        std::string dataName(this->volWithSol->getDataSet()->GetCellData()->GetArrayName(i));
        writer->writeSolutionField(dataName, nodeName, RealDouble,&cellData[0u]);
      }
    }
  } 
}

  //if (!partitions.empty())
  //{
  //  for (int i = 0; i < partitions.size(); ++i)
  //  {
  //    if (partitions[i])
  //    {
  //      delete partitions[i];
  //      partitions[i] = nullptr;
  //    }
  //  }
  //}
  //if (!virtualCellsOfPartitions.empty())
  //{
  //  for (int i = 0; i < virtualCellsOfPartitions.size(); ++i)
  //  {
  //    auto it = virtualCellsOfPartitions[i].begin();
  //    while (it != virtualCellsOfPartitions[i].end())
  //    {
  //      if (it->second)
  //      {
  //        delete it->second;
  //        it->second = nullptr;
  //      }
  //      ++it;
  //    }
  //  }
  //}
  //if (!surfacePartitions.empty())
  //{
  //  for (int i = 0; i < surfacePartitions.size(); ++i)
  //  {
  //    if (surfacePartitions[i])
  //    {
  //      delete surfacePartitions[i];
  //      surfacePartitions[i] = nullptr;
  //    }
  //  }
  //}
  //if (!virtualCellsOfSurfPartitions.empty())
  //{
  //  for (int i = 0; i < virtualCellsOfSurfPartitions.size(); ++i)
  //  {
  //    auto it = virtualCellsOfSurfPartitions[i].begin();
  //    while (it != virtualCellsOfSurfPartitions[i].end())
  //    {
  //      if (it->second)
  //      {
  //        delete it->second;
  //        it->second = nullptr;
  //      }
  //      ++it;
  //    }
  //  }
  //}
  //if (!patchesOfSurfacePartitions.empty())
  //{
  //  for (int i = 0; i < patchesOfSurfacePartitions.size(); ++i)
  //  {
  //    auto it = patchesOfSurfacePartitions[i].begin();
  //    while (it != patchesOfSurfacePartitions[i].end())
  //    {
  //      if (it->second)
  //      {
  //        delete it->second;
  //        it->second = nullptr;
  //      }
  //      ++it;
  //    }
  //  }
  //}
  //if (!virtualCellsOfPatchesOfSurfacePartitions.empty())
  //{
  //  for (int i = 0; i < virtualCellsOfPatchesOfSurfacePartitions.size(); ++i)
  //  {
  //    auto it = virtualCellsOfPatchesOfSurfacePartitions[i].begin();
  //    while (it != virtualCellsOfPatchesOfSurfacePartitions[i].end())
  //    {
  //      auto it1 = it->second.begin();
  //      while (it1 != it->second.end())
  //      {
  //        if (it1->second)
  //        {
  //          delete it1->second;
  //          it1->second = nullptr;
  //        }
  //        ++it1;
  //      }    
  //    ++it;
  //    }
  //  }
  //}
//  auto it = sharedPatchNodes[0].begin();
//  std::cout << "Partition 0 : " << std::endl;
//  while (it != sharedPatchNodes[0].end())
//  {
//    std::cout << "Patch " << it->first << " shared nodes with ";
//    auto it1 = it->second.begin();
//    while (it1 != it->second.end())
//    {
//      std::cout << "proc " << it1->first << ", patch ";
//      auto it2 = it1->second.begin();
//      while (it2 != it1->second.end())
//      {
//        std::cout << it2->first << " : \n";
//        for (int k = 0; k < it2->second.size(); ++k)
//          std::cout << it2->second[k] << " ";
//        std::cout << std::endl;
//        ++it2;
//      }
//      std::cout << std::endl;
//      ++it1;
//    }
//    std::cout << std::endl;  
//    ++it;
//  } 
//  std::cout << "pconn for proc 0, patch 3" << std::endl;
//  printVec(surfPconns[0][3]);
//  std::cout << std::endl;
  //for (int i = 0; i < partSurf->GetNumberOfPoints(); ++i)
  //{
  //  // get point from part surf
  //  double point[3];
  //  partSurf->GetPoint(i,point);
  //  // params to find cell with locator
  //  double closestPoint[3];
  //  vtkIdType id;
  //  int subid;
  //  double minDist2;
  //  // look for point in full surf cells
  //  fullSurfCellLocator->FindClosestPoint(point, closestPoint, genCell, id, subid, minDist2);  
  //  if (minDist2 > this->searchTolerance)
  //  {
  //    // find cells in part surf using point from part surf not found in full surf
  //    partSurf->GetPointCells(i, patchCellIds);
  //    for (int j = 0; j < patchCellIds->GetNumberOfIds(); ++j)
  //    {
  //      int cellToDelete = static_cast<int>(patchCellIds->GetId(j));
  //      cellsToDelete.insert(cellToDelete);
  //    }
  //    partSurf->DeletePoint(i);
  //    numToDelete += 1;
  //  }
  //} 
