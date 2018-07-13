#include <PreRocPrepDriver.H>
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

PreRocPrepDriver::PreRocPrepDriver(std::string& fname, int numPartitions, int numPatches)
{
  // load full volume mesh and create METIS partitions
  this->mesh = meshBase::Create(fname);
  this->partitions = meshBase::partition(mesh, numPartitions);
  // load stitched surf mesh with patch info
  this->stitchedSurf = meshBase::Create("stitchedSurf.vtu");
  this->AddGlobalCellIds(this->stitchedSurf);
  stitchedSurf->setContBool(0);

  std::vector<std::string> patchNoAndGlobalCellIds = {"patchNo","GlobalCellIds"};
  // create surfs partitions casted from vtp to vtu
  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
  this->surfacePartitions.resize(numPartitions);
  // initialize storage for virtual cells
  this->virtualCellsOfPartitions.resize(numPartitions);
  this->virtualCellsOfSurfPartitions.resize(numPartitions);
  // get all required local-global node maps and node/cell ids
  this->getGlobalIdsAndMaps(numPartitions,true);
  for (int i = 0; i < numPartitions; ++i)
  {
    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    appendFilter->AddInputData(
                    deleteInterPartitionSurface(
                      stitchedSurf, partitions[i]->extractSurface()));
    appendFilter->Update();
    unstructuredGrid = appendFilter->GetOutput();
    std::string basename("surfacePartition");
    basename += std::to_string(i);
    basename += ".vtu";
    // construct meshBase surface partition from vtkUnstructuredGrid
    meshBase* mbSurfPart = meshBase::Create(unstructuredGrid, basename);
    stitchedSurf->transfer(mbSurfPart, "Consistent Interpolation", patchNoAndGlobalCellIds, 1);
    //stitchedSurf->transfer(mbSurfPart,"Consistent Interpolation");
    mbSurfPart->write();
    this->surfacePartitions[i] = mbSurfPart;
    this->getGhostInformation(i,true);
    for (int j = 0; j < numPartitions; ++j)
    {
      this->getVirtualCells(i,j,true);
    }
  }
  this->getGlobalIdsAndMaps(numPartitions, false);
  for (int i = 0; i < numPartitions; ++i)
  {
    this->getGhostInformation(i,false);
    for (int j = 0; j < numPartitions; ++j)
    {
      this->getVirtualCells(i,j,false);
    }
  }  
  // create map of proc to patch to patch mesh
  this->extractPatches(numPatches);
}

void PreRocPrepDriver::AddGlobalCellIds(meshBase* _mesh)
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

void PreRocPrepDriver::getGlobalIdsAndMaps(int numPartitions, bool vol)
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

void PreRocPrepDriver::getVirtualCells(int me, int you, bool vol)
{
  if (vol && this->sharedNodes[me][you].size())
  {
    std::vector<int> virtuals(receivedCells[me][you].begin(), receivedCells[me][you].end());
    this->virtualCellsOfPartitions[me][you] = meshBase::extractSelectedCells(this->mesh,virtuals);
    std::stringstream ss;
    ss << "virtual" << "Vol" << "Of" << me << "from" << you << ".vtu";
    this->virtualCellsOfPartitions[me][you]->setFileName(ss.str());
    this->virtualCellsOfPartitions[me][you]->write();
  }
  else if (!vol && this->sharedSurfNodes[me][you].size())
  {
    std::vector<int> virtuals(receivedSurfCells[me][you].begin(), receivedSurfCells[me][you].end());
    this->virtualCellsOfSurfPartitions[me][you] 
      = meshBase::extractSelectedCells(this->stitchedSurf,virtuals);
    std::stringstream ss;
    ss << "virtual" << "Surf" << "Of" << me << "from" << you << ".vtu";
    this->virtualCellsOfSurfPartitions[me][you]->setFileName(ss.str());
    this->virtualCellsOfSurfPartitions[me][you]->write();
  }
}

void PreRocPrepDriver::getGhostInformation(int me, bool volOrSurf)
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

// TODO: need to get virtual info for patch and shared node stuff
//void PreRocPrepDriver::getSharedPatchInformation(int me)
//{
//  for (int i = 0; i < patchesOfSurfacePartitions.size(); ++i)
//  {
//
//  }
//}

void PreRocPrepDriver::getGhostInformation(int me, int you, bool hasShared, bool vol,
                                           vtkSmartPointer<vtkIdList> cellIdsList,
                                           vtkSmartPointer<vtkGenericCell> genCell)
{
  meshBase* meMesh = (vol ? partitions[me] : surfacePartitions[me]);

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

vtkSmartPointer<vtkPolyData>
PreRocPrepDriver::deleteInterPartitionSurface(meshBase* fullSurf, 
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
    double weights[3];
    double pcoords[3];
    double closestPoint[3];
    vtkIdType id;
    int subid;
    double minDist2;
    // look for point in full surf cells
    fullSurfCellLocator->FindClosestPoint(point, closestPoint, genCell, id, subid, minDist2);  
    if (minDist2 > 1e-9)//fullSurfCellLocator->FindCell(point, 0, genCell, pcoords, weights) == -1)
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
  vtkSmartPointer<vtkCleanPolyData> cleanPolyData 
    = vtkSmartPointer<vtkCleanPolyData>::New();
  cleanPolyData->SetInputData(partSurf);
  cleanPolyData->Update();
  return cleanPolyData->GetOutput();
}

void PreRocPrepDriver::extractPatches(int numPatches)
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
      meshBase* extractedPatch 
        = meshBase::extractSelectedCells(this->surfacePartitions[i], patchCellIds);
      std::stringstream ss;
      ss << "extractedPatch" << it->first << "FromProc" << i << ".vtu";
      extractedPatch->write(ss.str());
      patchesOfSurfacePartitions[i][it->first] = extractedPatch;  
      ++it;
    }
  } 
}


PreRocPrepDriver* PreRocPrepDriver::readJSON(json inputjson)
{
  std::string fname = inputjson["Grid To Partition"].as<std::string>();
  int numPartitions = inputjson["Number Of Partitions"].as<int>();
  int numPatches = inputjson["Number Of Patches"].as<int>();
  return new PreRocPrepDriver(fname, numPartitions, numPatches);
}

PreRocPrepDriver::~PreRocPrepDriver()
{
  if (mesh)
  {
    delete mesh;
    mesh = nullptr;
  }
  if (stitchedSurf)
  {
    delete stitchedSurf;
    stitchedSurf = nullptr;
  }
  if (!partitions.empty())
  {
    for (int i = 0; i < partitions.size(); ++i)
    {
      if (partitions[i])
      {
        delete partitions[i];
        partitions[i] = nullptr;
      }
    }
  }
  if (!virtualCellsOfPartitions.empty())
  {
    for (int i = 0; i < virtualCellsOfPartitions.size(); ++i)
    {
      auto it = virtualCellsOfPartitions[i].begin();
      while (it != virtualCellsOfPartitions[i].end())
      {
        if (it->second)
        {
          delete it->second;
          it->second = nullptr;
        }
        ++it;
      }
    }
  }
  if (!virtualCellsOfSurfPartitions.empty())
  {
    for (int i = 0; i < virtualCellsOfSurfPartitions.size(); ++i)
    {
      auto it = virtualCellsOfSurfPartitions[i].begin();
      while (it != virtualCellsOfSurfPartitions[i].end())
      {
        if (it->second)
        {
          delete it->second;
          it->second = nullptr;
        }
        ++it;
      }
    }
  }
  if (!surfacePartitions.empty())
  {
    for (int i = 0; i < surfacePartitions.size(); ++i)
    {
      if (surfacePartitions[i])
      {
        delete surfacePartitions[i];
        surfacePartitions[i] = nullptr;
      }
    }
  }
  if (!patchesOfSurfacePartitions.empty())
  {
    for (int i = 0; i < patchesOfSurfacePartitions.size(); ++i)
    {
      auto it = patchesOfSurfacePartitions[i].begin();
      while (it != patchesOfSurfacePartitions[i].end())
      {
        if (it->second)
        {
          delete it->second;
          it->second = nullptr;
        }
        ++it;
      }
    }
  }
}

  //// holds cell indices for each patch
  //std::vector<vtkSmartPointer<vtkIdTypeArray>> selectionIdsVec(numPatches);
  //for (int i = 0; i < numPatches; ++i)
  //{
  //  // creating selectionId container and push into vec
  //  vtkSmartPointer<vtkIdTypeArray> selectionIds = vtkSmartPointer<vtkIdTypeArray>::New();
  //  selectionIds->SetNumberOfComponents(1);
  //  selectionIdsVec[i] = selectionIds;
  //}
  //// getting patch number array from stitched surf
  //vtkSmartPointer<vtkDataArray> patchNumbers 
  //  = this->stitchedSurf->getDataSet()->GetCellData()->GetArray("patchNo");
  //for (int i = 0; i < patchNumbers->GetNumberOfTuples(); ++i)
  //{
  //  double tmp[1];
  //  patchNumbers->GetTuple(i,tmp);
  //  int patchNo = (int) tmp[0];
  //  selectionIdsVec[patchNo-1]->InsertNextValue(i);
  //}
  //// size vector of surface partition patches
  //this->surfacePartitionPatches.resize(numPatches);
  //for (int i = 0; i < numPatches; ++i)
  //{
  //  
  //  meshBase* extractedSurf 
  //    = meshBase::extractSelectedCells(this->stitchedSurf->getDataSet(), selectionIdsVec[i]);
  //  std::stringstream ss;
  //  ss << "extractedPatch" << i << ".vtu";
  //  extractedSurf->write(ss.str());
  //  this->surfacePartitionPatches[i] = extractedSurf;
  //}



      //std::cout << "Num nodes proc " << i << " shares with " << "proc " << j;
      //std::cout << " : " << this->sharedNodes[i][j].size() 
      //          << " " << this->sharedNodes[j][i].size() << std::endl;
      //std::cout << "Num nodes proc " << i << " sends to " << "proc " << j;
      //std::cout << " : " << this->sentNodes[i][j].size() << std::endl;  
      //std::cout << "Num nodes proc " << j << " receives from " << "proc " << i;
      //std::cout << " : " << this->receivedNodes[j][i].size() << std::endl;  
      //
      //std::vector<int> virtuals(receivedCells[i][j].size());
      //auto it = receivedCells[i][j].begin();
      //int c=0;
      //while (it != receivedCells[i][j].end())
      //{
      //  virtuals[c] = *it;
      //  ++c;
      //  ++it;
      //}
      //meshBase* tmp = meshBase::extractSelectedCells(mesh,virtuals);
      //std::stringstream ss;
      //ss << "virtualOf" << i << "from" << j << ".vtu";
      //tmp->write(ss.str());
      //delete tmp;
    
