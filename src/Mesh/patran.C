#include <meshBase.H>
#include <patran.H>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <fstream>
#include <algorithm>
#include <vtkCellData.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkGenericCell.h>


using namespace PATRAN;

//////////////////////////////////
// patran class 
//////////////////////////////////


// todo are patran indices 0-indexed or 1-indexed?

// Write title card
void patran::write25(std::ofstream& outputStream)
{
  std::cout << "writing packet 25..." << std::endl;

  // writer header
  std::string defaultHeader = "25       0       0       1       0       0       0       0       0";
  outputStream << defaultHeader;
  outputStream << std::endl;

  // write title
  outputStream << outFnameNeu;
  outputStream << std::endl;
}

// Write summary data
void patran::write26(std::ofstream& outputStream)
{
  std::cout << "writing packet 26..." << std::endl;
  
  // write header
  outputStream << "26";
  for (int i = 0; i < 8; i++)
  {
    switch(i)
    {
      case 2:
        outputStream << std::setw(8) << std::right << "1";
        break;
      case 3:
        outputStream << std::setw(8) << std::right << volMeshBase->getDataSet()->GetNumberOfPoints();
        break;
      case 4:
        outputStream << std::setw(8) << std::right << volMeshBase->getDataSet()->GetNumberOfCells();
        break;
      default:
        outputStream << std::setw(8) << std::right << "0";
        break;
    }
  }
  outputStream << std::endl;

  // Write time stamp information
  std::time_t t = std::time(0);   // get time now
  std::tm* now = std::localtime(&t);
  outputStream << std::setw(12) << std::left << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday;
  outputStream << std::endl;  
}

// Write node data
void patran::write1(std::ofstream& outputStream)
{
  std::cout << "writing packet 1..." << std::endl;

  for (int iNode = 0; iNode < volMeshBase->getDataSet()->GetNumberOfPoints(); iNode++)
  {
    // Header card
    // packet no
    std::string defaultHeader = " 1";
    outputStream << defaultHeader;
    // node id
    outputStream << std::setw(8) << std::right << iNode+1;
    // IV, default to 0
    outputStream << std::setw(8) << std::right << "0";
    // KC = 2
    outputStream << std::setw(8) << std::right << "2";
    // fill in zeros
    for (int i = 0; i < 5; i++)
    {
      outputStream << std::setw(8) << std::right << "0";
    }
    outputStream << std::endl;

    // Data Card 1
    // Get node coordinates
    double coords[3];
    volMeshBase->getDataSet()->GetPoint(iNode, coords);
    for (int i = 0; i < 3; i++)
    {
      char buffer[17];
      snprintf(buffer,17,"%16.9E",coords[i]);
      outputStream << buffer;
      //outputStream << std::setw(16) << std::left << printf("%14.12E",coords[i]);
      //outputStream << std::setw(16) << std::right << std::scientific << std::setprecision(9) << coords[i];
    }
    outputStream << std::endl;

    // Data Card 2
    // currently defaults to ICF=1, GTYPE=G for first two
    defaultHeader = "1G";
    outputStream << defaultHeader;
    // degrees of freedom
    outputStream << std::setw(8) << std::right << "6";
    // node configuration
    outputStream << std::setw(8) << std::right << "0";
    // CID (coordinate frame for analysis results)
    outputStream << std::setw(8) << std::right << "0";
    // Write two white space
    outputStream << std::setw(2) << std::right << "";
    // PSPC, permament single point constraing flags
    outputStream << std::setw(6) << std::right << "000000";

    outputStream << std::endl;
  }
}

// Write element data
void patran::write2(std::ofstream& outputStream)
{

  std::cout << "writing packet 2..." << std::endl;
  
  for (int iCell = 0; iCell < volMeshBase->getDataSet()->GetNumberOfCells(); iCell++)
  {
    // ONLY WRITE FOR TETRAHEDRAL ELEMENTS:
    double elemType[1];
    volMeshBase->getDataSet()->GetCellData()->GetArray("elemType")
                                                  ->GetTuple(iCell, elemType);
    if (elemType[0] == VTK_TETRA)
    {
      // Header card
      // packet no
      std::string defaultHeader = " 2";
      outputStream << defaultHeader;
      // cell id
      outputStream << std::setw(8) << std::right << iCell+1;
      // IV, default to TETRAHEDRAL = 5
      outputStream << std::setw(8) << std::right << "5";
      // KC = 1 + (NODES + 9)/10 + (N1 + 4)/5
      // Default values for TETRAHEDRAL:
      // NODES = 4, N1 = 0 ---> 1+(4+9)/10+(0+4)/5 = 3.1 ---> floor(3.1) = 3??
      // Wrote 2 below as in example files...
      outputStream << std::setw(8) << std::right << "2";
      // fill in zeros (defaulting N1 and N2 to zero)
      for (int i = 0; i < 5; i++)
      {
        outputStream << std::setw(8) << std::right << "0";
      }
      outputStream << std::endl;
    
      // Data Card 1
      // Default to 4 nodes per TETRAHEDRAL
      outputStream << std::setw(8) << std::right << "4";
      // CONFIG (element type), set it to TETRAHEDRAL (5)
      // It is used in Packet 4, which Rocstar doesn't not support, so
      // probably unneeded.
      outputStream << std::setw(8) << std::right << "0";
      // PID, Property ID (unused)
      outputStream << std::setw(8) << std::right << "0";
      // CEID, Congruent Element ID (unused)
      outputStream << std::setw(8) << std::right << "0";
      // Material Orientation Angles (unused)
      for (int i = 0; i < 3; i++)
      {
        char buffer[17];
        snprintf(buffer,17,"%16.9E",0.0);
        outputStream << buffer;
        //outputStream << std::setw(16) << std::right << std::scientific << std::setprecision(9) << 0.0;
      }
      outputStream << std::endl;
  
      // Data Card 2
      // Get tetrahedral cell points
      // Pts beyond 4th can be omitted (don't need to write buffer zeros)
      vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
      volMeshBase->getDataSet()->GetCellPoints(iCell, cellPtIds);
      for (int iPt = 0; iPt < cellPtIds->GetNumberOfIds(); iPt++)
      {
        outputStream << std::setw(8) << std::right << cellPtIds->GetId(iPt)+1;
      }
      outputStream << std::endl;
  
      // Data Card 3 is ignored
    }
  }
}

// Write distributed loads
void patran::write6(std::ofstream& outputStream)
{
  std::cout << "writing packet 6..." << std::endl;
  
  vtkSmartPointer<vtkIdList> facePtIds;
  vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New(); 
  vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New();
  
  // building cell locator for looking up patch number in remeshed surface mesh
  vtkSmartPointer<vtkCellLocator> surfCellLocator = surfMeshBase->buildLocator(); 
  // maximum number of vertices per face (to be found in proceeding loop)
  int nVerticesPerFaceMax = 0;
  // maximum number of faces per cell (to be found in proceeding loop)
  int nFacesPerCellMax = 0; 

  // Accumulate map of patch numbers of each boundary node
  //std::map<int,std::vector<int>> boundaryNodeId2PatchNo;

  std::map<std::vector<int>, std::pair<int,int>, sortIntVec_compare> faceMap;
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
      vtkCell* faceObj = genCell1->GetFace(j);
      bool shared = 0;
      int numVerts = faceObj->GetNumberOfPoints();
      nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts ? numVerts : nVerticesPerFaceMax);
      facePtIds = faceObj->GetPointIds(); 
      volMeshBase->getDataSet()->GetCellNeighbors(i, facePtIds, sharedCellPtIds); 
      std::vector<int> facePntIds(numVerts);
      for (int k = 0; k < numVerts; ++k)
      {
        facePntIds[k] = faceObj->GetPointId(k)+1;
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
        faceObj->GetPoints()->GetPoint(0,p1);
        faceObj->GetPoints()->GetPoint(1,p2);
        faceObj->GetPoints()->GetPoint(2,p3);
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

        // Assign patch numbers to all nodes in the face
        for (int iFaceNode = 0; iFaceNode < facePtIds->GetNumberOfIds(); iFaceNode++)
        {
          boundaryNodeId2PatchNo[facePtIds->GetId(iFaceNode)].push_back(patchNo[0]);
        }

        // Header card
        // packet number
        outputStream << std::setw(2) << std::right << "6";
        // element number
        outputStream << std::setw(8) << std::right << i;
        // default load set = 1
        outputStream << std::setw(8) << std::right << "1";
        // defalt KC = 2
        outputStream << std::setw(8) << std::right << "2";
        // Fill in remaining zeros
        for (int i = 0; i < 5; i++)
        {
          outputStream << std::setw(8) << std::right << "0";
        }
        outputStream << std::endl;

        // Data Card 1
        // default face value
        outputStream << std::setw(3) << std::right << "110"; // unused by Rocfrac
        
        // load direction
        //outputStream << std::setw(6) << std::right << distLoadDirMap[patchNo[0]];
        outputStream << std::setw(6) << std::right << "000000"; // ICOMP, unused by Rocfrac
        
        // face nodes
        outputStream << std::setw(8) << std::right << face2nodes[j];
        
        // face number
        //outputStream << std::setw(2) << std::right << j+1; // 1-indexed in PATRAN
        outputStream << std::setw(2) << std::right << "1"; // NFE, unused by Rocfrac
        outputStream << std::endl;

        // Data Card 2
        // card number
        char buffer[17];
        //std::cout << "patchNo[0] = " << patchNo[0] << std::endl;
        //std::cout << "faceTypeMap[patchNo[0]] = " << faceTypeMap[patchNo[0]] << std::endl;
        snprintf(buffer,17,"%16.9E",(double) faceTypeMap[patchNo[0]]);
        outputStream << buffer;

        //outputStream << std::setw(16) << std::setprecision(9) << std::scientific << (double) faceTypeMap[patchNo[0]];
        outputStream << std::endl;
      }
    }
  }
}

//// Write node forces
//void patran::write7(std::ofstream& outputStream)
//{
//
//  std::cout << "writing packet 7..." << std::endl;
//
//  for (int iNode = 0; iNode < volMeshBase->getDataSet()->GetNumberOfPoints(); iNode++)
//  {
//    std::cout << "inode = " << iNode << std::endl;
//    vtkSmartPointer<vtkIdList> pointCellIds = vtkSmartPointer<vtkIdList>::New();
//    // Todo how to set nodes that belong to two patches?
//    volMeshBase->getDataSet()->GetPointCells(iNode, pointCellIds);
//
//    for (int jCell = 0; jCell < pointCellIds->GetNumberOfIds(); jCell++)
//    {
//      vtkIdType cellId = pointCellIds->GetId(jCell);
//      vtkSmartPointer<vtkIdList> facePtIds;
//      vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
//      vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New(); 
//      vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New(); 
//      std::map<std::vector<int>, std::pair<int,int>, sortIntVec_compare> faceMap;
//      // building cell locator for looking up patch number in remeshed surface mesh
//      vtkSmartPointer<vtkCellLocator> surfCellLocator = surfMeshBase->buildLocator(); 
//      // maximum number of vertices per face (to be found in proceeding loop)
//      int nVerticesPerFaceMax = 0;
//      // maximum number of faces per cell (to be found in proceeding loop)
//      int nFacesPerCellMax = 0; 
//
//      volMeshBase->getDataSet()->GetCell(cellId, genCell1);
//      for (int kFace = 0; kFace < genCell1->GetNumberOfFaces(); ++kFace)
//      {
//        vtkCell* faceObj = genCell1->GetFace(kFace);
//        bool shared = 0;
//        int numVerts = faceObj->GetNumberOfPoints();
//        nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts ? numVerts : nVerticesPerFaceMax);
//        facePtIds = faceObj->GetPointIds(); 
//        volMeshBase->getDataSet()->GetCellNeighbors(jCell, facePtIds, sharedCellPtIds);
//
//        if (!(sharedCellPtIds->GetNumberOfIds()))
//        {
//          double p1[3], p2[3], p3[3];
//          faceObj->GetPoints()->GetPoint(0,p1);
//          faceObj->GetPoints()->GetPoint(1,p2);
//          faceObj->GetPoints()->GetPoint(2,p3);
//          double faceCenter[3];
//          for (int k = 0; k < numVerts; ++k)
//          {
//            faceCenter[k] = (p1[k] + p2[k] + p3[k])/3.0;
//          } 
//          vtkIdType closestCellId;
//          int subId;
//          double minDist2;
//          double closestPoint[3];
//          // find closest point and closest cell to faceCenter
//          surfCellLocator->FindClosestPoint(
//            faceCenter, closestPoint, genCell2,closestCellId,subId,minDist2);
//          double patchNo[1];
//          surfMeshBase->getDataSet()->GetCellData()->GetArray("patchNo")
//                                                    ->GetTuple(closestCellId, patchNo);
//
//
//
//          // Header card
//          // packet number
//          outputStream << std::setw(2) << std::left << "7";
//          // node number
//          outputStream << std::setw(8) << std::left << iNode;
//          // default load set = 1
//          outputStream << std::setw(8) << std::left << "1";
//          // default KC = 2
//          outputStream << std::setw(8) << std::left << "2";
//          // Fill in remaining zeros
//          for (int i = 0; i < 5; i++)
//          {
//            outputStream << std::setw(8) << std::left << "0";
//          }
//  
//          // Data Card 1
//          // coordinate frame id, set default to 0
//          outputStream << std::setw(8) << std::left << "0";  // default CID, unused by Rocfrac
//          //outputStream << std::setw(6) << std::left << nodeForceDirMap[patchNo[0]];
//          outputStream << std::setw(6) << std::left << "100000" // default ICOMP, unused by Rocfrac
//  
//          // Data Card 2
//          // card number
//          outputStream << std::setw(16) << std::setprecision(9) << std::scientific << (double)nodeForceCardMap[patchNo[0]];
//          outputStream << std::endl;
//        }
//      }
//    }
//  }
//}

bool patran::comparePatch(int i, int j)
{
  auto posi = find(nppVec.begin(), nppVec.end(), i);
  auto posj = find(nppVec.begin(), nppVec.end(), j);
  //std::cout << "i = " << i << std::endl;
  //std::cout << "posi = " << posi-nppVec.begin() << std::endl;
  //std::cout << "j = " << j << std::endl;
  //std::cout << "posj = " << posj-nppVec.begin() << std::endl;
  if (posi-nppVec.begin() < posj-nppVec.begin())
  {
    return true;
  }
  else
  {
    return false;
  }
}

// Write node displacements
void patran::write8(std::ofstream& outputStream)
{

  std::cout << "writing packet 8..." << std::endl;

  for (auto nodeItr = boundaryNodeId2PatchNo.begin(); nodeItr != boundaryNodeId2PatchNo.end(); ++nodeItr)
//  for (int iNode = 0; iNode < volMeshBase->getDataSet()->GetNumberOfPoints(); iNode++)
  {
    int iNode = nodeItr->first;

    vtkSmartPointer<vtkIdList> pointCellIds = vtkSmartPointer<vtkIdList>::New();
    // Todo how to set nodes that belong to two patches?
    volMeshBase->getDataSet()->GetPointCells(iNode, pointCellIds);

    // Header card
    // packet number
    outputStream << std::setw(2) << std::right << "8";
    // node number
    outputStream << std::setw(8) << std::right << iNode+1;
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
    //std::cout << "original order = " << std::endl;
    for (auto patchItr = boundaryNodeId2PatchNo[iNode].begin(); patchItr != boundaryNodeId2PatchNo[iNode].end(); patchItr++)
    {
      //std::cout << *patchItr << std::endl;
      patchNoVec.push_back(*patchItr);
    }
    // sort patches
    int outFlag = 0;
    int patchVal = -1;
    for (auto patchItr = patchNoVec.begin(); patchItr != patchNoVec.end()-1; ++patchItr)
    {
      int firstVal = *patchItr;
      if (patchVal == -1)
      {
        patchVal = *patchItr;
      }
      else
      {
        //if (*patchItr != patchVal)
        //{
        //  std::cout << "comparing new ... " << *patchItr << " with old " << patchVal << std::endl;
        //}
        if (comparePatch(*patchItr,patchVal))
        {
          //std::cout << "----------------" << std::endl;
          //std::cout << "old Patch = " << patchVal << std::endl;
          //std::cout << "new patch = " << *patchItr << std::endl;
          patchVal = *patchItr;
        }
      }
      //int secondVal = *(patchItr+1);
      //if (firstVal != secondVal)
      //{
      //  outFlag = 1;
      //  //std::cout << "comparing " << firstVal << " and " << secondVal << std::endl;
      //  if (!comparePatch(firstVal,secondVal))
      //  {
      //   //std::cout << "swapping..." << std::endl;
      //   iter_swap(patchItr,patchItr+1);
      //  }
      //}
    }
    //std::sort(patchNoVec.begin(),patchNoVec.end(),comparePatch);
    // put in array
    int patchNo[boundaryNodeId2PatchNo[iNode].size()];
    //int iPatch = 0;
    //if (outFlag)
    //{
    //  std::cout << "new order = ";
    //}
    //for (auto patchItr = patchNoVec.begin(); patchItr != patchNoVec.end(); ++patchItr)
    //{
    //  if (outFlag)
    //  {
    //    std::cout << *patchItr << " ";
    //  }
    //  patchNo[iPatch] = *patchItr;
    //  iPatch++;
    //}
    //if (outFlag)
    //{
    //  std::cout << std::endl;
    //}
//

    // Assign patches according to node patch preference ordering
    // sort patch numbers by preference
    //std::vector<int> patchNoTmp;
    //iPatch = 0;
    //for (auto nppItr = nppVec.begin(); nppItr != nppVec.end(); ++nppItr)
    //{
    //  if (patchNoTmp.size() == 0)
    //  {
    //    patchNoTmp.push_back(*nppItr);
    //  }
    //  else
    //  {
    //    for (auto pnItr = patchNoTmp.begin(); pnItr != patchNoTmp.end(); ++pnItr)
    //    {
    //      auto it = std::find(pnnVec.begin(), pnnVec.end(), *pnItr);
    //      if (it)
    //      posCurrent = ;
    //      posInsert;
    //    }
    //  }
    //}

    patchNo[0] = patchVal;

    // Data Card 1
    // coordinate frame id, set default to 0
    outputStream << std::setw(8) << std::right << "0";  // default CID, unused by Rocfrac
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

    outputStream << std::setw(6) << std::right << structuralFlag + thermalFlag + meshMotionFlag + "000";
    //outputStream << std::setw(6) << std::left << "100000" // default ICOMP, unused by Rocfrac
    outputStream << std::endl;

    // Data Card 2
    // card number
    for (int i = 0; i < nFlags; i++)
    {
      char buffer[17];
      snprintf(buffer,17,"%16.9E",(double) nodeTypeMap[patchNo[0]]);
      outputStream << buffer;

      //outputStream << std::setw(16) << std::setprecision(9) << std::scientific << (double) nodeTypeMap[patchNo[0]];
    }
    outputStream << std::endl;
  }
}

//    for (int jCell = 0; jCell < pointCellIds->GetNumberOfIds(); jCell++)
//    {
//      vtkIdType cellId = pointCellIds->GetId(jCell);
//      vtkSmartPointer<vtkIdList> facePtIds;
//      vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
//      vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New(); 
//      vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New(); 
//      std::map<std::vector<int>, std::pair<int,int>, sortIntVec_compare> faceMap;
//      // building cell locator for looking up patch number in remeshed surface mesh
//      vtkSmartPointer<vtkCellLocator> surfCellLocator = surfMeshBase->buildLocator(); 
//      // maximum number of vertices per face (to be found in proceeding loop)
//      int nVerticesPerFaceMax = 0;
//      // maximum number of faces per cell (to be found in proceeding loop)
//      int nFacesPerCellMax = 0; 
//
//      volMeshBase->getDataSet()->GetCell(cellId, genCell1);
//      for (int kFace = 0; kFace < genCell1->GetNumberOfFaces(); ++kFace)
//      {
//        vtkCell* faceObj = genCell1->GetFace(kFace);
//        bool shared = 0;
//        int numVerts = faceObj->GetNumberOfPoints();
//        nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts ? numVerts : nVerticesPerFaceMax);
//        facePtIds = faceObj->GetPointIds(); 
//        volMeshBase->getDataSet()->GetCellNeighbors(jCell, facePtIds, sharedCellPtIds);
//
//        if (!(sharedCellPtIds->GetNumberOfIds()))
//        {
//          double p1[3], p2[3], p3[3];
//          faceObj->GetPoints()->GetPoint(0,p1);
//          faceObj->GetPoints()->GetPoint(1,p2);
//          faceObj->GetPoints()->GetPoint(2,p3);
//          double faceCenter[3];
//          for (int k = 0; k < numVerts; ++k)
//          {
//            faceCenter[k] = (p1[k] + p2[k] + p3[k])/3.0;
//          } 
//          vtkIdType closestCellId;
//          int subId;
//          double minDist2;
//          double closestPoint[3];
//          // find closest point and closest cell to faceCenter
//          surfCellLocator->FindClosestPoint(
//            faceCenter, closestPoint, genCell2,closestCellId,subId,minDist2);
//          double patchNo[1];
//          surfMeshBase->getDataSet()->GetCellData()->GetArray("patchNo")
//                                                    ->GetTuple(closestCellId, patchNo);
//
//
//
//          // Header card
//          // packet number
//          outputStream << std::setw(2) << std::left << "8";
//          // node number
//          outputStream << std::setw(8) << std::left << iNode;
//          // default load set = 1
//          outputStream << std::setw(8) << std::left << "1";
//          // default KC = 2
//          outputStream << std::setw(8) << std::left << "2";
//          // Fill in remaining zeros
//          for (int i = 0; i < 5; i++)
//          {
//            outputStream << std::setw(8) << std::left << "0";
//          }
//  
//          // Data Card 1
//          // coordinate frame id, set default to 0
//          outputStream << std::setw(8) << std::left << "0";  // default CID, unused by Rocfrac
//          //outputStream << std::setw(6) << std::left << nodeDispDirMap[patchNo[0]];
//          outputStream << std::setw(6) << std::left << "100000" // default ICOMP, unused by Rocfrac
//  
//          // Data Card 2
//          // card number
//          outputStream << std::setw(16) << std::setprecision(9) << std::scientific << nodeDispCardMap[patchNo[0]];
//          outputStream << std::endl;
//        }
//      }
//    }
//  }
//}

//// Write node tempratures
//void patran::write10(std::ofstream& outputStream)
//{
//
//  std::cout << "writing packet 10..." << std::endl;
//
//  // IS THIS EVEN SUPPORTED BY ROCFRAC??
//}

// Write element data
void patran::write99(std::ofstream& outputStream)
{
  outputStream << "99       0       0       1       0       0       0       0       0";
  outputStream << std::endl;
}

// constructor from meshBase object
patran::patran(const std::shared_ptr<meshBase> _fullMesh, const std::string _inFnameVtk,
  const std::string _outFnameNeu, std::map<int, int> _faceTypeMap,
  std::map<int, int> _nodeTypeMap, std::map<int, bool> _nodeStructuralMap,
  std::map<int, bool> _nodeMeshMotionMap, std::map<int, bool> _nodeThermalMap,
  std::vector<int> _nppVec)
{

  std::cout << "writing patran file..." << std::endl;

  fullMesh = _fullMesh;
  inFnameVtk = _inFnameVtk;
  outFnameNeu = _outFnameNeu;
  faceTypeMap = _faceTypeMap;
  nodeTypeMap = _nodeTypeMap;
  nodeStructuralMap = _nodeStructuralMap;
  nodeMeshMotionMap = _nodeMeshMotionMap;
  nodeThermalMap = _nodeThermalMap;
  nppVec = _nppVec;

  // Todo is this the right ordering?
  //face2nodes[0] = "01110000";
  //face2nodes[1] = "10110000";
  //face2nodes[2] = "11010000";
  //face2nodes[3] = "11100000";

  //face2nodes[0] = "11100000";
  //face2nodes[1] = "11010000";
  //face2nodes[2] = "01110000";
  //face2nodes[3] = "10110000";

  face2nodes[0] = "11100000";
  face2nodes[1] = "11010000";
  face2nodes[2] = "10110000";
  face2nodes[3] = "01110000";

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

  // instantiate output stream for writing
  std::ofstream outputStream(outFnameNeu);
  if(!outputStream.good()) 
  {
    std::cout << "Cannot open file " << outFnameNeu << std::endl;
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
  volMeshBase->write();  // write to file

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
  surfMeshBase->write();  // write to file

  // write title card
  write25(outputStream);
  // write summary data
  write26(outputStream);
  // write node data
  write1(outputStream);
  // write element data
  write2(outputStream);
  // write distributed loads
  write6(outputStream);
  // write all node BCs
  write8(outputStream);
  // write end of file flag
  write99(outputStream);





  // write node temperatures
  //write10(outputStream);
  // write node forces
  //write7(outputStream);
}