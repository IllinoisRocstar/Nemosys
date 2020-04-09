#include <algorithm>

#include "patran.H"

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkGenericCell.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

//////////////////////////////////
// patran class 
//////////////////////////////////


// Write title card (packet 25)
void PATRAN::patran::write25(std::ofstream &outputStream) const
{
  std::cout << "writing packet 25..." << std::endl;

  // write header card
  std::string defaultHeader = "25       0       0       1       0       0       0       0       0";
  outputStream << defaultHeader;
  outputStream << std::endl;

  // write title card
  outputStream << outFnameNeu;
  outputStream << std::endl;
}

// Write summary data (packet 26)
void PATRAN::patran::write26(std::ofstream &outputStream) const
{
  std::cout << "writing packet 26..." << std::endl;

  // write header card
  outputStream << "26";
  for (int i = 0; i < 8; i++)
  {
    switch (i)
    {
      case 2:
        outputStream << std::setw(8) << std::right << "1";
        break;
      case 3:
        outputStream << std::setw(8) << std::right
                     << volMeshBase->getDataSet()->GetNumberOfPoints();
        break;
      case 4:
        outputStream << std::setw(8) << std::right
                     << volMeshBase->getDataSet()->GetNumberOfCells();
        break;
      default:
        outputStream << std::setw(8) << std::right << "0";
        break;
    }
  }
  outputStream << std::endl;

  // Write time stamp info card
  std::time_t t = std::time(nullptr);
  std::tm *now = std::localtime(&t);
  outputStream << std::setw(12) << std::left << (now->tm_year + 1900) << '-'
               << (now->tm_mon + 1) << '-' << now->tm_mday;
  outputStream << std::endl;
}

// Write node data (packet 1)
void PATRAN::patran::write1(std::ofstream &outputStream) const
{
  std::cout << "writing packet 1..." << std::endl;

  for (int iNode = 0;
       iNode < volMeshBase->getDataSet()->GetNumberOfPoints(); iNode++)
  {
    // Write header card
    // packet no
    outputStream << " 1";
    // node id (1-indexed)
    outputStream << std::setw(8) << std::right << iNode + 1;
    // IV, default to 0
    outputStream << std::setw(8) << std::right << "0";
    // KC, default to 2
    outputStream << std::setw(8) << std::right << "2";
    // fill in zeros for remaining
    for (int i = 0; i < 5; i++)
    {
      outputStream << std::setw(8) << std::right << "0";
    }
    outputStream << std::endl;

    // Write data card 1
    // Get node coordinates
    double coords[3];
    volMeshBase->getDataSet()->GetPoint(iNode, coords);
    for (double coord : coords)
    {
      char buffer[17];
      snprintf(buffer, 17, "%16.9E", coord);
      outputStream << buffer;
    }
    outputStream << std::endl;

    // Write data card 2
    // currently defaults to ICF=1, GTYPE=G for first two
    outputStream << "1G";
    // degrees of freedom
    outputStream << std::setw(8) << std::right << "6";
    // node configuration
    outputStream << std::setw(8) << std::right << "0";
    // CID (coordinate frame for analysis results)
    outputStream << std::setw(8) << std::right << "0";
    // Write two white spaces
    outputStream << std::setw(2) << std::right << "";
    // PSPC, permanent single point constraining flags
    outputStream << std::setw(6) << std::right << "000000";
    outputStream << std::endl;
  }
}

// Write element data (packet 2)
void PATRAN::patran::write2(std::ofstream &outputStream) const
{
  std::cout << "writing packet 2..." << std::endl;

  for (int iCell = 0;
       iCell < volMeshBase->getDataSet()->GetNumberOfCells(); iCell++)
  {
    // Note: only tetrahedral elements are written

    // Get element type
    double elemType[1];
    volMeshBase->getDataSet()->GetCellData()->GetArray("elemType")
        ->GetTuple(iCell, elemType);
    if (elemType[0] == VTK_TETRA)
    {
      // Write header card
      // packet no
      std::string defaultHeader = " 2";
      outputStream << defaultHeader;
      // cell id
      outputStream << std::setw(8) << std::right << iCell + 1;
      // IV, default to 5 (tetrahedral)
      outputStream << std::setw(8) << std::right << "5";
      // KC, default to 2
      outputStream << std::setw(8) << std::right << "2";
      // fill in zeros (defaulting N1 and N2 to zero)
      for (int i = 0; i < 5; i++)
      {
        outputStream << std::setw(8) << std::right << "0";
      }
      outputStream << std::endl;

      // Write data card 1
      // Default to 4 nodes per tetrahedral
      outputStream << std::setw(8) << std::right << "4";
      // CONFIG, element type (unused)
      outputStream << std::setw(8) << std::right << "0";
      // PID, Property ID (unused)
      outputStream << std::setw(8) << std::right << "0";
      // CEID, Congruent Element ID (unused)
      outputStream << std::setw(8) << std::right << "0";
      // Material Orientation Angles (unused)
      for (int i = 0; i < 3; i++)
      {
        char buffer[17];
        snprintf(buffer, 17, "%16.9E", 0.0);
        outputStream << buffer;
      }
      outputStream << std::endl;

      // Write data Card 2
      // Get tetrahedral cell points
      // Pts beyond 4th can be omitted (don't need to write buffer zeros)
      vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
      volMeshBase->getDataSet()->GetCellPoints(iCell, cellPtIds);
      for (int iPt = 0; iPt < cellPtIds->GetNumberOfIds(); iPt++)
      {
        outputStream << std::setw(8) << std::right << cellPtIds->GetId(iPt) + 1;
      }
      outputStream << std::endl;

      // Data Card 3 is ignored
    }
  }
}

// Write distributed loads (packet 6)
void PATRAN::patran::write6(std::ofstream &outputStream)
{
  std::cout << "writing packet 6..." << std::endl;

  vtkSmartPointer<vtkIdList> facePtIds;
  vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New();

  // building cell locator for looking up patch number in remeshed surface mesh
  vtkSmartPointer<vtkCellLocator> surfCellLocator = surfMeshBase->buildLocator();
  // maximum number of vertices per face (to be found in proceeding loop)
  vtkIdType nVerticesPerFaceMax = 0;
  // maximum number of faces per cell (to be found in proceeding loop)
  int nFacesPerCellMax = 0;

  std::map<std::vector<nemId_t>,
           std::pair<nemId_t, nemId_t>,
           sortNemId_tVec_compare> faceMap;
  for (nemId_t i = 0; i < volMeshBase->getNumberOfCells(); ++i)
  {
    // get cell i 
    volMeshBase->getDataSet()->GetCell(i, genCell1);
    // get faces, find cells sharing it. if no cell shares it, 
    // use the locator of the surfMeshBase to find the patch number
    int numFaces = genCell1->GetNumberOfFaces();
    nFacesPerCellMax = (nFacesPerCellMax < numFaces
                        ? numFaces
                        : nFacesPerCellMax);
    for (int j = 0; j < numFaces; ++j)
    {
      vtkCell *faceObj = genCell1->GetFace(j);
      // bool shared = false;
      vtkIdType numVerts = faceObj->GetNumberOfPoints();
      nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts
                             ? numVerts
                             : nVerticesPerFaceMax);
      facePtIds = faceObj->GetPointIds();

      volMeshBase->getDataSet()->GetCellNeighbors(i, facePtIds,
                                                  sharedCellPtIds);
      std::vector<nemId_t> facePntIds(numVerts);
      for (vtkIdType k = 0; k < numVerts; ++k)
      {
        facePntIds[k] = faceObj->GetPointId(k) + 1;
      }
      if (sharedCellPtIds->GetNumberOfIds())
      {
        faceMap.insert(
            std::pair<std::vector<nemId_t>, std::pair<nemId_t, nemId_t>>(
                facePntIds,
                std::make_pair(
                    i + 1,
                    sharedCellPtIds->GetId(0) + 1
                )
            )
        );
      }
      else
      {
        double p1[3], p2[3], p3[3];
        faceObj->GetPoints()->GetPoint(0, p1);
        faceObj->GetPoints()->GetPoint(1, p2);
        faceObj->GetPoints()->GetPoint(2, p3);
        double faceCenter[3];
        for (nemId_t k = 0; k < numVerts; ++k)
        {
          faceCenter[k] = (p1[k] + p2[k] + p3[k]) / 3.0;
        }
        vtkIdType closestCellId;
        int subId;
        double minDist2;
        double closestPoint[3];
        // find closest point and closest cell to faceCenter
        surfCellLocator->FindClosestPoint(
            faceCenter, closestPoint, genCell2, closestCellId, subId, minDist2);
        double patchNo[1];
        surfMeshBase->getDataSet()->GetCellData()->GetArray("patchNo")
            ->GetTuple(closestCellId, patchNo);
        faceMap.insert(
            std::pair<std::vector<nemId_t>, std::pair<nemId_t, nemId_t>>(
                facePntIds,
                std::make_pair(
                    i + 1,
                    -1 * patchNo[0]
                )
            )
        );

        // Assign patch numbers to all nodes in the face
        for (vtkIdType iFaceNode = 0;
             iFaceNode < facePtIds->GetNumberOfIds(); iFaceNode++)
        {
          boundaryNodeId2PatchNo[facePtIds->GetId(iFaceNode)].push_back(
              static_cast<int>(patchNo[0]));
        }

        // Write header card
        // packet number
        outputStream << std::setw(2) << std::right << "6";
        // element number, 1-indexed
        outputStream << std::setw(8) << std::right << i + 1;
        // default load set = 1
        outputStream << std::setw(8) << std::right << "1";
        // default KC = 2
        outputStream << std::setw(8) << std::right << "2";
        // Fill in remaining zeros
        for (int i = 0; i < 5; i++)
        {
          outputStream << std::setw(8) << std::right << "0";
        }
        outputStream << std::endl;

        // Write data card 1
        // default face value
        outputStream << std::setw(3) << std::right
                     << "110"; // unused by Rocfrac

        // load direction
        outputStream << std::setw(6) << std::right
                     << "000000"; // ICOMP, unused by Rocfrac

        // face nodes
        outputStream << std::setw(8) << std::right << face2nodes[j];

        // face number
        outputStream << std::setw(2) << std::right
                     << "1"; // NFE, unused by Rocfrac
        outputStream << std::endl;

        // Write data card 2
        // card number
        char buffer[17];
        snprintf(buffer, 17, "%16.9E", (double) faceTypeMap[static_cast<int>(patchNo[0])]);
        outputStream << buffer;

        outputStream << std::endl;
      }
    }
  }
}

// compare patches in nppVec to determine preference
bool PATRAN::patran::comparePatch(int i, int j)
{
  auto posi = find(nppVec.begin(), nppVec.end(), i);
  auto posj = find(nppVec.begin(), nppVec.end(), j);
  return posi - nppVec.begin() < posj - nppVec.begin();
}

// Write node BCs: structural, mesh motion, and thermal (packet 8)
void PATRAN::patran::write8(std::ofstream &outputStream)
{
  std::cout << "writing packet 8..." << std::endl;

  for (const auto &nodeItr : boundaryNodeId2PatchNo)
  {
    nemId_t iNode = nodeItr.first;

    vtkSmartPointer<vtkIdList> pointCellIds = vtkSmartPointer<vtkIdList>::New();
    volMeshBase->getDataSet()->GetPointCells(iNode, pointCellIds);

    // Write header card
    // packet number
    outputStream << std::setw(2) << std::right << "8";
    // node number, 1-indexed
    outputStream << std::setw(8) << std::right << iNode + 1;
    // default load set = 1
    outputStream << std::setw(8) << std::right << "1";
    // default KC = 2
    outputStream << std::setw(8) << std::right << "2";
    // Fill in remaining zeros
    for (int i = 0; i < 5; i++)
    {
      outputStream << std::setw(8) << std::right << "0";
    }
    outputStream << std::endl;

    // Get patch no
    std::vector<int> patchNoVec;
    for (int &patchItr : boundaryNodeId2PatchNo[iNode])
      patchNoVec.push_back(patchItr);

    // sort patches
    // int outFlag = 0;
    int patchVal = -1;
    for (auto patchItr = patchNoVec.begin();
         patchItr != patchNoVec.end() - 1; ++patchItr)
    {
      // int firstVal = *patchItr;
      if (patchVal == -1)
        patchVal = *patchItr;
      else
        if (comparePatch(*patchItr, patchVal))
          patchVal = *patchItr;
    }
    int patchNo[1] = {patchVal};

    // Write data card 1
    // coordinate frame id, set default to 0
    outputStream << std::setw(8) << std::right
                 << "0";  // default CID, unused by Rocfrac
    std::string structuralFlag, thermalFlag, meshMotionFlag;
    int nFlags = 0;
    if (nodeStructuralMap[patchNo[0]])
    {
      structuralFlag = "1";
      nFlags++;
    }
    else
    {
      structuralFlag = "0";
    }
    if (nodeThermalMap[patchNo[0]])
    {
      thermalFlag = "1";
      nFlags++;
    }
    else
    {
      thermalFlag = "0";
    }
    if (nodeMeshMotionMap[patchNo[0]])
    {
      meshMotionFlag = "1";
      nFlags++;
    }
    else
    {
      meshMotionFlag = "0";
    }

    outputStream << std::setw(6) << std::right
                 << structuralFlag << thermalFlag << meshMotionFlag << "000";
    outputStream << std::endl;

    // Write data card 2
    // card number
    for (int i = 0; i < nFlags; i++)
    {
      char buffer[17];
      snprintf(buffer, 17, "%16.9E", (double) nodeTypeMap[patchNo[0]]);
      outputStream << buffer;
    }
    outputStream << std::endl;
  }
}

// Write "end of file" line (packet 99)
void PATRAN::patran::write99(std::ofstream &outputStream) const
{
  outputStream
      << "99       0       0       1       0       0       0       0       0";
  outputStream << std::endl;
}

// constructor from meshBase object
PATRAN::patran::patran(const std::shared_ptr<meshBase> _fullMesh,
                       const std::string &_inFnameVtk,
                       const std::string &_outFnameNeu,
                       const std::map<int, int> &_faceTypeMap,
                       const std::map<int, int> &_nodeTypeMap,
                       const std::map<int, bool> &_nodeStructuralMap,
                       const std::map<int, bool> &_nodeMeshMotionMap,
                       const std::map<int, bool> &_nodeThermalMap,
                       const std::vector<int> &_nppVec)
{
  std::cout << "Writing Patran file in Rocfrac input format..." << std::endl;

  fullMesh = _fullMesh;
  inFnameVtk = _inFnameVtk;
  outFnameNeu = _outFnameNeu;
  faceTypeMap = _faceTypeMap;
  nodeTypeMap = _nodeTypeMap;
  nodeStructuralMap = _nodeStructuralMap;
  nodeMeshMotionMap = _nodeMeshMotionMap;
  nodeThermalMap = _nodeThermalMap;
  nppVec = _nppVec;

  // Face/Node ordering convention: face2Nodes[faceNo] = nodes used; 1 = used, 0 = unused
  face2nodes[0] = "11010000";
  face2nodes[1] = "01110000";
  face2nodes[2] = "10110000";
  face2nodes[3] = "11100000";

  //if (inFnameVtk.find(".vt") != -1) 
  //{
  //  // get trimmed filename
  //  size_t lastindex = inFnameVtk.find_last_of("."); 
  //}
  //else
  //{
  //  std::cout << "File must be in VTK format." << std::endl;
  //  exit(1);
  //}

  // instantiate output stream for writing
  std::ofstream outputStream(outFnameNeu);
  if (!outputStream.good())
  {
    std::cout << "Cannot open file " << outFnameNeu << std::endl;
    exit(1);
  }

  // declare array to store element type
  vtkSmartPointer<vtkDataArray> elementTypeArray = vtkSmartPointer<vtkIdTypeArray>::New();
  elementTypeArray->SetNumberOfComponents(1);
  elementTypeArray->SetName("elemType");

  // get element types
  double tmp[1];  // storing inserted tuples
  for (vtkIdType iCell = 0;
       iCell < fullMesh->getDataSet()->GetNumberOfCells(); iCell++)
  {
    if (fullMesh->getDataSet()->GetCellType(iCell) != VTK_TETRA
        && fullMesh->getDataSet()->GetCellType(iCell) != VTK_TRIANGLE)
    {
      std::cout << "Error: only triangle and tetrahedral elements supported."
                << std::endl;
      exit(1);
    }
    tmp[0] = fullMesh->getDataSet()->GetCellType(iCell);
    elementTypeArray->InsertNextTuple(tmp);
  }
  fullMesh->getDataSet()->GetCellData()->AddArray(elementTypeArray);

  // threshold to extract only volume cells
  vtkSmartPointer<vtkThreshold> volThreshold = vtkSmartPointer<vtkThreshold>::New();
  volThreshold->SetInputData(fullMesh->getDataSet());
  volThreshold->SetInputArrayToProcess(0, 0, 0,
                                       vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                       "elemType");
  volThreshold->ThresholdBetween(VTK_TETRA, VTK_TETRA);
  volThreshold->Update();
  // store output in unstructured grid
  vtkSmartPointer<vtkUnstructuredGrid> volUG = volThreshold->GetOutput();
  // create meshbase object
  volMeshBase = meshBase::CreateShared(volUG, "extractedVolume.vtu");
  //volMeshBase->write();  // write to file

  // threshold to extract only surface cells
  vtkSmartPointer<vtkThreshold> surfThreshold = vtkSmartPointer<vtkThreshold>::New();
  surfThreshold->SetInputData(fullMesh->getDataSet());
  surfThreshold->SetInputArrayToProcess(0, 0, 0,
                                        vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                        "elemType");
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
  //surfMeshBase->write();  // write to file

  // Write Patran file packets
  std::cout << "Writing Patran file..." << std::endl;

  // write title card
  write25(outputStream);

  // write summary data
  write26(outputStream);

  // write node data
  write1(outputStream);

  // write element data
  write2(outputStream);

  // write distributed load BCs
  write6(outputStream);

  // write all node BCs
  write8(outputStream);

  // write end of file flag
  write99(outputStream);
}
