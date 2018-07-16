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

RocPartCommGenDriver::RocPartCommGenDriver(std::shared_ptr<meshBase> _mesh, 
                                           std::shared_ptr<meshBase> _remeshedSurf, 
                                           std::shared_ptr<meshBase> _volWithSol,
                                           std::shared_ptr<meshBase> _surfWithSol,
                                           int numPartitions)
{
  std::cout << "RocPartCommGenDriver created\n";
  this->mesh = _mesh;
  // load stitched surf mesh with patch info
  this->remeshedSurf = _remeshedSurf;
  this->volWithSol = _volWithSol;
  this->surfWithSol = _surfWithSol; 
  this->execute(numPartitions);
}


RocPartCommGenDriver::RocPartCommGenDriver(std::string& volname, std::string& surfname, int numPartitions)
{
  std::cout << "RocPartCommGenDriver created\n";
  // load full volume mesh and create METIS partitions
  this->mesh = meshBase::CreateShared(volname);
  this->remeshedSurf = meshBase::CreateShared(surfname);
  this->execute(numPartitions);
}

void RocPartCommGenDriver::execute(int numPartitions)
{
  remeshedSurf->setContBool(0);
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
    //remeshedSurf->transfer(mbSurfPart,"Consistent Interpolation");
    //mbSurfPart->write();
    //this->surfacePartitions[i] = mbSurfPart;
    this->surfacePartitions[i]->write();
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
      // get virtual cells of each volume partition (t4:real)
      this->getVirtualCells(i,j,true);
    }
    // transfer data from original volume with solution if it's provided
    if (this->volWithSol)
    {
      this->volWithSol->setContBool(0);
      this->volWithSol->transfer(this->partitions[i].get(),"Consistent Interpolation");
    }
    // write cgns for vol partition
    this->writeCgns(i,1);
  }
  // get global ids and maps for surface partitions
  this->getGlobalIdsAndMaps(numPartitions, false);
  for (int i = 0; i < numPartitions; ++i)
  {
    // get ghost information for each surface partition
    this->getGhostInformation(i,false);
    for (int j = 0; j < numPartitions; ++j)
    {
      // get virtual cells for each surface partition
      this->getVirtualCells(i,j,false);
    }
  }  
  // create map of proc to patch to patch mesh
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

void RocPartCommGenDriver::getVirtualCells(int me, int you, bool vol)
{
  if (vol && this->sharedNodes[me][you].size())
  {
    std::vector<int> virtuals(receivedCells[me][you].begin(), receivedCells[me][you].end());
    this->virtualCellsOfPartitions[me][you] = 
      meshBase::CreateShared(meshBase::extractSelectedCells(this->mesh.get(),virtuals));
    std::stringstream ss;
    ss << "virtual" << "Vol" << "Of" << me << "from" << you << ".vtu";
    this->virtualCellsOfPartitions[me][you]->setFileName(ss.str());
    this->virtualCellsOfPartitions[me][you]->write();
  }
  else if (!vol && this->sharedSurfNodes[me][you].size())
  {
    std::vector<int> virtuals(receivedSurfCells[me][you].begin(), receivedSurfCells[me][you].end());
    this->virtualCellsOfSurfPartitions[me][you] 
      = meshBase::CreateShared(meshBase::extractSelectedCells(this->remeshedSurf.get(),virtuals));
    std::stringstream ss;
    ss << "virtual" << "Surf" << "Of" << me << "from" << you << ".vtu";
    this->virtualCellsOfSurfPartitions[me][you]->setFileName(ss.str());
    this->virtualCellsOfSurfPartitions[me][you]->write();
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
      // get the cell in me for point extraction
      meMesh->getDataSet()->GetCell(localCellId, genCell);
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
        std::cout << patchProcGlobNodeIdsMaps[it->first][me].size() << std::endl;
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
  vtkSmartPointer<vtkIdList> patchCellIds = vtkSmartPointer<vtkIdList>::New();
  std::unordered_set<int> cellsToDelete;
  int numToDelete = 0;
  for (int i = 0; i < partSurf->GetNumberOfPoints(); ++i)
  {
    // get point from part surf
    double point[3];
    partSurf->GetPoint(i,point);
    // params to find cell with locator
    double closestPoint[3];
    vtkIdType id;
    int subid;
    double minDist2;
    // look for point in full surf cells
    fullSurfCellLocator->FindClosestPoint(point, closestPoint, genCell, id, subid, minDist2);  
    if (minDist2 > 1e-9)
    {
      // find cells in part surf using point from part surf not found in full surf
      partSurf->GetPointCells(i, patchCellIds);
      for (int j = 0; j < patchCellIds->GetNumberOfIds(); ++j)
      {
        int cellToDelete = static_cast<int>(patchCellIds->GetId(j));
        cellsToDelete.insert(cellToDelete);
      }
      partSurf->DeletePoint(i);
      numToDelete += 1;
    }
  } 
  std::unordered_set<int>::iterator it = cellsToDelete.begin();
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
      std::vector<int> patchCellIds(it->second.begin(), it->second.end());  
      this->patchesOfSurfacePartitions[i][it->first] 
        = meshBase::CreateShared(
            meshBase::extractSelectedCells(this->surfacePartitions[i].get(), patchCellIds));
      std::stringstream ss;
      ss << "extractedPatch" << it->first << "OfProc" << i << ".vtu";
      this->patchesOfSurfacePartitions[i][it->first]->setFileName(ss.str());//write(ss.str());
      this->patchesOfSurfacePartitions[i][it->first]->write();
      std::cout << ss.str() << std::endl;
      ++it;
    }
  }
  this->virtualCellsOfPatchesOfSurfacePartitions.resize(surfacePartitions.size());
  for (int i = 0; i < this->virtualCellsOfSurfPartitions.size(); ++i)
  {
    std::map<int, std::vector<int>> virtualPatchPartitionCellMap;
    auto it = virtualCellsOfSurfPartitions[i].begin();
    while (it != virtualCellsOfSurfPartitions[i].end())
    {
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
        std::vector<int> virtualPatchCellIds(it1->second.begin(), it1->second.end());
        this->virtualCellsOfPatchesOfSurfacePartitions[i][it->first][it1->first]
          = meshBase::CreateShared(
              meshBase::extractSelectedCells(it->second.get(), virtualPatchCellIds));
        std::stringstream ss;
        ss << "extractedVirtualPatch" << it1->first 
          << "Of" << i << "FromProc" << it->first  << ".vtu";
        this->virtualCellsOfPatchesOfSurfacePartitions[i][it->first][it1->first]->setFileName(ss.str());
        this->virtualCellsOfPatchesOfSurfacePartitions[i][it->first][it1->first]->write();
        std::cout << ss.str() << std::endl;
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

void RocPartCommGenDriver::writeCgns(int proc, int type)
{
  std::stringstream ss;
  ss << "fluid_00.000000_000" << proc << ".cgns";
  std::unique_ptr<cgnsWriter> writer 
    = std::unique_ptr<cgnsWriter>(new cgnsWriter(ss.str(), "fluid", 3, 3));
  // define elementary information 
  writer->setUnits(Kilogram,Meter,Second,Kelvin,Degree);
  // baseitrname, nTstep, timeVal 
  writer->setBaseItrData("TimeIterValues", 1, 0.0);
  // setting zone iter data
  writer->setZoneItrData("ZoneIterativeData", "GridCoordinatesPointers", "FlowSolutionPointers");
  ss.str("");
  ss.clear();
  ss << 0 << proc+1 << 0 << type;
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
  // write solution data
  //if (this->volWithSol)
  //{
  //  // begin with point data
  //  for (int i = 0; i < this->volWithSol->getDataSet()->GetPointData()->GetNumberOfArrays())
  //  {
  //    std::vector<double> pointData;
  //    this->volWithSol->getPointDataArray

  //  }

  //  

  //} 
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
