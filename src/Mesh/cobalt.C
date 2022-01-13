/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#include "Mesh/cobalt.H"

// VTK
#include <vtkCellTypes.h>
#include <vtkCellData.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkGenericCell.h>
#include <vtkUnstructuredGrid.h>

#include "AuxiliaryFunctions.H"

using nemAux::operator-;  // for vector subtraction.

//////////////////////////////////
// cobalt class 
//////////////////////////////////

// constructor from meshBase object
COBALT::cobalt::cobalt(const std::shared_ptr<meshBase> fullMesh,
                       const std::string &_inFnameVtk,
                       const std::string &_outFnameCgr,
                       const std::string &_outFnameCgi)
{
  inFnameVtk = _inFnameVtk;     // vtk input file
  outFnameCgr = _outFnameCgr;   // cobalt output grid file
  outFnameCgi = _outFnameCgi;   // cobalt output patch mapping file

  // declare array to store element type
  vtkSmartPointer<vtkDataArray> elementTypeArray
      = vtkSmartPointer<vtkIdTypeArray>::New();
  elementTypeArray->SetNumberOfComponents(1);
  elementTypeArray->SetName("elemType");

  // get element types
  double tmp[1];  // storing inserted tuples
  for (int iCell = 0;
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
  // create meshBase object
  volMeshBase = meshBase::CreateShared(volUG, "extractedVolume.vtu");
  //volUG_mb->write();  // write to file

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
  // geometry filter to polyData object
  vtkSmartPointer<vtkGeometryFilter> geomFilter = vtkGeometryFilter::New();
  geomFilter->SetInputData(surfUG);
  geomFilter->Update();
  vtkSmartPointer<vtkPolyData> surfPD = geomFilter->GetOutput();
  // create meshBase object
  surfMeshBase = meshBase::CreateShared(surfPD, "extractedSurface.vtp");
  //surfPD_mb->write();  // write to file
}

// writes COBALT mesh data into file
void COBALT::cobalt::write() const
{
  std::ofstream outputStream(outFnameCgr);
  if (!outputStream.good())
  {
    std::cout << "Cannot open file " << outFnameCgr << std::endl;
    exit(1);
  }

  if (!surfMeshBase)
  {
    std::cout << "surface mesh is empty!" << std::endl;
    exit(1);
  }
  auto surfPatchNoArr =
      surfMeshBase->getDataSet()->GetCellData()->GetArray("patchNo");
  if (!surfPatchNoArr) {
    surfPatchNoArr =
        surfMeshBase->getDataSet()->GetCellData()->GetArray("PhysGrpId");
    if (!surfPatchNoArr) {
      std::cout << "surface mesh must have cell array with name \"patchNo\" or "
                   "\"PhysGrpId\" denoting patch number."
                << std::endl;
      exit(1);
    }
  }

  vtkSmartPointer<vtkIdList> facePtIds;
  vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New();
  // Signed to represent patch numbers as negative integers
  std::map<std::vector<nemId_t>,
           std::pair<nemId_t, long long int>,
           sortNemId_tVec_compare> faceMap;
  // building cell locator for looking up patch number in remeshed surface mesh
  vtkSmartPointer<vtkStaticCellLocator> surfCellLocator = surfMeshBase->buildStaticCellLocator();
  // maximum number of vertices per face (to be found in proceeding loop)
  vtkIdType nVerticesPerFaceMax = 0;
  // maximum number of faces per cell (to be found in proceeding loop)
  int nFacesPerCellMax = 0;

  for (vtkIdType i = 0; i < volMeshBase->getNumberOfCells(); ++i)
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
      vtkCell *face = genCell1->GetFace(j);
      // bool shared = false;
      vtkIdType numVerts = face->GetNumberOfPoints();
      nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts
                             ? numVerts
                             : nVerticesPerFaceMax);
      facePtIds = face->GetPointIds();
      volMeshBase->getDataSet()->GetCellNeighbors(i, facePtIds,
                                                  sharedCellPtIds);
      std::vector<nemId_t> facePntIds(numVerts);
      for (vtkIdType k = 0; k < numVerts; ++k)
      {
        facePntIds[k] = face->GetPointId(k) + 1;
      }
      if (sharedCellPtIds->GetNumberOfIds())
      {
        faceMap.emplace(facePntIds, std::pair<nemId_t, long long int>(
                                        i + 1, sharedCellPtIds->GetId(0) + 1));
      }
      else
      {
        double bounds[6];
        face->GetBounds(bounds);
        auto cellIds = vtkSmartPointer<vtkIdList>::New();
        surfCellLocator->FindCellsWithinBounds(bounds, cellIds);
        vtkIdType closestCellId;
        if (cellIds->GetNumberOfIds() == 0) {
          std::cerr << "No surface cell found." << std::endl;
          continue;
        } else {
          closestCellId = cellIds->GetId(0);
          if (cellIds->GetNumberOfIds() > 1) {
            double p1[3], p2[3], p3[3];
            face->GetPoints()->GetPoint(0, p1);
            face->GetPoints()->GetPoint(1, p2);
            face->GetPoints()->GetPoint(2, p3);
            std::vector<double> faceCenter(3);
            for (vtkIdType k = 0; k < numVerts; ++k) {
              faceCenter[k] = (p1[k] + p2[k] + p3[k]) / 3.0;
            }
            double minDist = nemAux::l2_Norm(
                faceCenter - surfMeshBase->getCellCenter(cellIds->GetId(0)));
            for (vtkIdType k = 1; k < cellIds->GetNumberOfIds(); ++k) {
              auto candidateCellCenter =
                  surfMeshBase->getCellCenter(cellIds->GetId(k));
              auto dist = nemAux::l2_Norm(faceCenter - candidateCellCenter);
              if (dist < minDist) {
                minDist = dist;
                closestCellId = cellIds->GetId(k);
              }
            }
          }
        }
        double patchNo = surfPatchNoArr->GetComponent(closestCellId, 0);
        faceMap.emplace(facePntIds,
                        std::pair<nemId_t, long long int>(i + 1, -1 * patchNo));
      }
    }
  }

  std::map<nemId_t, nemId_t> patchMap;
  for (nemId_t i = 0; i < surfMeshBase->getNumberOfCells(); ++i)
  {
    double patchNo[1];
    surfPatchNoArr->GetTuple(i, patchNo);
    patchMap.insert(std::pair<nemId_t, nemId_t>(patchNo[0], i));
  }

  // write patch mapping file
  writePatchMap(outFnameCgi, patchMap);
  // write cobalt file
  outputStream << 3 << "   " << 1 << "  " << patchMap.size() << std::endl;
  outputStream << volMeshBase->getNumberOfPoints() << " " << faceMap.size()
               << " " << volMeshBase->getNumberOfCells() << " "
               << nVerticesPerFaceMax << " " << nFacesPerCellMax << std::endl;
  for (nemId_t i = 0; i < volMeshBase->getNumberOfPoints(); ++i)
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


void COBALT::cobalt::writePatchMap(const std::string &mapFile,
                                   const std::map<nemId_t, nemId_t> &patchMap) const
{
  std::ofstream outputStream(mapFile);
  if (!outputStream.good())
  {
    std::cout << "Error opening file " << mapFile << std::endl;
    exit(1);
  }
  writePatchMap(outputStream, patchMap);
}


void COBALT::cobalt::writePatchMap(std::ofstream &outputStream,
                                   const std::map<nemId_t, nemId_t> &patchMap) const
{
  outputStream << patchMap.size() << std::endl;
  outputStream << patchMap.size() << std::endl;
  auto it = patchMap.begin();
  size_t normPatchNo = 1;
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
