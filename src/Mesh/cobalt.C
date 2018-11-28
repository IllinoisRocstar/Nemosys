#include <meshBase.H>
#include <cobalt.H>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <fstream>
#include <algorithm>
#include <vtkCellData.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkGenericCell.h>


using namespace COBALT;

//////////////////////////////////
// cobalt class 
//////////////////////////////////

// constructor from meshBase object
cobalt::cobalt(const std::shared_ptr<meshBase> fullMesh, const std::string _inFnameVtk, const std::string _outFnameCgr, const std::string _outFnameCgi)
{
  inFnameVtk = _inFnameVtk;
  outFnameCgr = _outFnameCgr;
  outFnameCgi = _outFnameCgi;
  if (inFnameVtk.find(".vt") != -1) 
  {
    // get trimmed filename
    size_t lastindex = inFnameVtk.find_last_of("."); 
  }
  else
  {
    std::cout << "File must be in VTK format." << std::endl;
    exit(1);
  }

  // declare array to store element type
  vtkSmartPointer<vtkDataArray> elementTypeArray = vtkSmartPointer<vtkIdTypeArray>::New();;
  elementTypeArray->SetNumberOfComponents(1);
  elementTypeArray->SetName("elemType");

  // get element types
  double tmp[1];  // storing inserted tuples
  for (int iCell = 0; iCell < fullMesh->getDataSet()->GetNumberOfCells(); iCell++)
  {
    if (fullMesh->getDataSet()->GetCellType(iCell) != VTK_TETRA &&
        fullMesh->getDataSet()->GetCellType(iCell) != VTK_TRIANGLE)
    {
      std::cout << "Only triangle and tetrahedral elements supported." << std::endl;
      exit(1);
    }
    tmp[0] = fullMesh->getDataSet()->GetCellType(iCell);
    elementTypeArray->InsertNextTuple(tmp);
  }
  fullMesh->getDataSet()->GetCellData()->AddArray(elementTypeArray);

  // threshold to get only volume cells
  vtkSmartPointer<vtkThreshold> volThreshold = vtkSmartPointer<vtkThreshold>::New();
  volThreshold->SetInputData(fullMesh->getDataSet());
  volThreshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "elemType");
  volThreshold->ThresholdBetween(VTK_TETRA, VTK_TETRA);
  volThreshold->Update();
  // store output in unstructured grid
  vtkSmartPointer<vtkUnstructuredGrid> volUG = volThreshold->GetOutput();
  // create meshbase object
  volMeshBase = meshBase::CreateShared(volUG, "extractedVolume.vtu");
  //volUG_mb->write();  // write to file

  // threshold to get only surface cells
  vtkSmartPointer<vtkThreshold> surfThreshold = vtkSmartPointer<vtkThreshold>::New();
  surfThreshold->SetInputData(fullMesh->getDataSet());
  surfThreshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "elemType");
  surfThreshold->ThresholdBetween(VTK_TRIANGLE, VTK_TRIANGLE);
  surfThreshold->Update();
  // store output in unstructured grid
  vtkSmartPointer<vtkUnstructuredGrid> surfUG = surfThreshold->GetOutput();
  // geometry filter to polydata object
  vtkSmartPointer<vtkGeometryFilter> geomFilter = vtkGeometryFilter::New();
  geomFilter->SetInputData(surfUG);
  geomFilter->Update();
  vtkSmartPointer<vtkPolyData> surfPD = geomFilter->GetOutput();
  // create meshbase object
  surfMeshBase = meshBase::CreateShared(surfPD, "extractedSurface.vtp");
  //surfPD_mb->write();  // write to file
}

// writes COBALT mesh data into file
void cobalt::write()
{

  std::ofstream outputStream(outFnameCgr);
  if(!outputStream.good()) 
  {
    std::cout << "Cannot open file " << outFnameCgr << std::endl;
    exit(1);
  }

  if (!surfMeshBase) 
  {
    std::cout << "surface mesh is empty!" << std::endl;
    exit(1);
  }
  if (surfMeshBase->IsArrayName("patchNo", 1) == -1)
  {
    std::cout << "surface mesh must have patchNo cell array" << std::endl;
    exit(1);
  } 
  vtkSmartPointer<vtkIdList> facePtIds;
  vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New(); 
  vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New(); 
  std::map<std::vector<int>, std::pair<int,int>, sortIntVec_compare> faceMap;
  // building cell locator for looking up patch number in remeshed surface mesh
  vtkSmartPointer<vtkCellLocator> surfCellLocator = surfMeshBase->buildLocator(); 
  // maximum number of vertices per face (to be found in proceeding loop)
  int nVerticesPerFaceMax = 0;
  // maximum number of faces per cell (to be found in proceeding loop)
  int nFacesPerCellMax = 0; 

  for (int i = 0; i < volMeshBase->getNumberOfCells(); ++i)
  {
    // get cell i 
    volMeshBase->getDataSet()->GetCell(i, genCell1);
    // get faces, find cells sharing it. if no cell shares it, 
    // use the locator of the surfMeshBase to find the patch number
    int numFaces = genCell1->GetNumberOfFaces();
    nFacesPerCellMax = (nFacesPerCellMax < numFaces ? numFaces : nFacesPerCellMax);
    for (int j = 0; j < numFaces; ++j)
    {
      vtkCell* face = genCell1->GetFace(j);
      bool shared = 0;
      int numVerts = face->GetNumberOfPoints();
      nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts ? numVerts : nVerticesPerFaceMax);
      facePtIds = face->GetPointIds(); 
      volMeshBase->getDataSet()->GetCellNeighbors(i, facePtIds, sharedCellPtIds); 
      std::vector<int> facePntIds(numVerts);
      for (int k = 0; k < numVerts; ++k)
      {
        facePntIds[k] = face->GetPointId(k)+1;
      }
      //std::cout << "sharedCellPtIds->GetNumberOfIds() = " << sharedCellPtIds->GetNumberOfIds() << std::endl;
      if (sharedCellPtIds->GetNumberOfIds())
      {
        faceMap.insert(std::pair<std::vector<int>, std::pair<int,int>>
                        (facePntIds, std::make_pair(i+1, (int) sharedCellPtIds->GetId(0)+1))); 
      }
      else
      {
        double p1[3], p2[3], p3[3];
        face->GetPoints()->GetPoint(0,p1);
        face->GetPoints()->GetPoint(1,p2);
        face->GetPoints()->GetPoint(2,p3);
        double faceCenter[3];
        for (int k = 0; k < numVerts; ++k)
        {
          faceCenter[k] = (p1[k] + p2[k] + p3[k])/3.0;
        } 
        vtkIdType closestCellId;
        int subId;
        double minDist2;
        double closestPoint[3];
        // find closest point and closest cell to faceCenter
        surfCellLocator->FindClosestPoint(
          faceCenter, closestPoint, genCell2,closestCellId,subId,minDist2);
        double patchNo[1];
        surfMeshBase->getDataSet()->GetCellData()->GetArray("patchNo")
                                                  ->GetTuple(closestCellId, patchNo);
        faceMap.insert(std::pair<std::vector<int>, std::pair<int,int>>      
                  (facePntIds, std::make_pair(i+1, (int) -1*patchNo[0])));
      }
    }
  }
  
  std::map<int,int> patchMap;
  for (int i = 0; i < surfMeshBase->getNumberOfCells(); ++i)
  {
    double patchNo[1];
    surfMeshBase->getDataSet()->GetCellData()->GetArray("patchNo")
                                              ->GetTuple(i, patchNo);
    patchMap.insert(std::pair<int,int>(patchNo[0],i));
  }

  // write patch mapping file
  writePatchMap(outFnameCgi, patchMap);
  // write cobalt file
  outputStream << 3 << "   " << 1 << "  " << patchMap.size() << std::endl;
  outputStream << volMeshBase->getNumberOfPoints() << " " << faceMap.size()
               << " " << volMeshBase->getNumberOfCells() << " "
               << nVerticesPerFaceMax << " " << nFacesPerCellMax << std::endl;
  for (int i = 0; i < volMeshBase->getNumberOfPoints(); ++i)
  {
    std::vector<double> pnt(volMeshBase->getPoint(i));
    outputStream << std::setw(21) << std::fixed << std::setprecision(15)
                 << pnt[0] << "   " << pnt[1] << "   " << pnt[2] << std::endl;
  }
  
  auto it = faceMap.begin();
  while (it != faceMap.end())
  {
    outputStream << it->first.size() << " ";
    auto faceIdIter = it->first.begin();
    while (faceIdIter != it->first.end())
    {
      outputStream << *faceIdIter << " ";
      ++faceIdIter;
    }
    outputStream << it->second.first << " " << it->second.second << std::endl;
    ++it;
  }

}

void cobalt::writePatchMap(const std::string& mapFile, const std::map<int,int>& patchMap)
{
  std::ofstream outputStream(mapFile);
  if (!outputStream.good())
  {
    std::cout << "Error opening file " << mapFile << std::endl;
    exit(1);
  }
  writePatchMap(outputStream, patchMap);
}

void cobalt::writePatchMap(std::ofstream& outputStream, const std::map<int,int>& patchMap)
{
  outputStream << patchMap.size() << std::endl;
  outputStream << patchMap.size() << std::endl;
  auto it = patchMap.begin();
  int normPatchNo = 1;
  while (it != patchMap.end())
  {
    for (int i = 0; i < 2; ++i)
    {
      outputStream << std::setw(2) << std::left << it->first << " ";
    }
    outputStream << std::setw(2) << std::left << normPatchNo << " ";
    outputStream << std::endl;
    ++it;
    normPatchNo++;
  }
}