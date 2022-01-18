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
#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/geoMeshBase.H"

#include <iostream>
#include <set>
#include <utility>

#ifdef HAVE_GMSH
#include <gmsh.h>
#endif

#include <vtkAppendFilter.h>

namespace NEM {
namespace MSH {

GmshInterface::GmshInterface() {
#ifdef HAVE_GMSH
  gmsh::initialize();
  gmsh::option::setNumber("General.Terminal", 1.0);  // Gmsh errors to stderr
  gmsh::option::setNumber("Mesh.SaveAll", 1);
  std::cout << "Gmsh initialized" << std::endl;
#endif
}

GmshInterface::~GmshInterface() {
#ifdef HAVE_GMSH
  gmsh::finalize();
  std::cout << "Gmsh finalized" << std::endl;
#endif
}

std::shared_ptr<GmshInterface> GmshInterface::instance = nullptr;
int GmshInterface::count = 0;

void GmshInterface::Initialize() {
  ++count;
  if (!instance) instance = std::shared_ptr<GmshInterface>(new GmshInterface());
}

void GmshInterface::Finalize() {
  --count;
  if (count == 0) instance.reset();
}

geoMeshBase::geoMeshBase()
    : geoMeshBase({vtkSmartPointer<vtkUnstructuredGrid>::New(), "", "", {}}) {}

geoMeshBase::geoMeshBase(GeoMesh inGeoMesh)
    : _geoMesh(std::move(inGeoMesh)), _angleThreshold(30.0 * M_PI / 180.0) {
  InitializeObjectBase();
#ifdef HAVE_GMSH
  GmshInterface::Initialize();
#endif
  std::cout << "geoMeshBase constructed" << std::endl;
}

geoMeshBase::~geoMeshBase() {
  this->ReferenceCount--;
#ifdef HAVE_GMSH
  GmshInterface::Finalize();
#endif
  std::cout << "geoMeshBase destructed" << std::endl;
}

void geoMeshBase::report(std::ostream &out) const {
  out << this->GetClassName() << ":\n";
  auto mesh = getGeoMesh().mesh;
  out << "Nodes:\t" << mesh->GetNumberOfPoints() << '\n';
  out << "Elements:\t" << mesh->GetNumberOfCells() << '\n';
}

void geoMeshBase::takeGeoMesh(geoMeshBase *otherGeoMesh) {
  _geoMesh = std::move(otherGeoMesh->_geoMesh);
  otherGeoMesh->_geoMesh = {
      vtkSmartPointer<vtkUnstructuredGrid>::New(), {}, {}, {}};
  otherGeoMesh->resetNative();
  resetNative();
}

// Only merges vtkDataSets. Need to add merge for sideSets, and geometry
void geoMeshBase::mergeGeoMesh(geoMeshBase *otherGeoMesh) {
  vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
  appendFilter->AddInputData(_geoMesh.mesh);
  appendFilter->AddInputData(otherGeoMesh->_geoMesh.mesh);
  appendFilter->MergePointsOn();
  appendFilter->Update();
  _geoMesh.mesh = appendFilter->GetOutput();
  otherGeoMesh->_geoMesh = {
      vtkSmartPointer<vtkUnstructuredGrid>::New(), {}, {}, {}};
  otherGeoMesh->resetNative();
}

// Helper functions for geoMeshBase::get{Point, Cell, Field}DataArrayCopy
namespace {

vtkSmartPointer<vtkAbstractArray> copyOrCreateCopy(
    vtkAbstractArray *orig, vtkAbstractArray *copy = nullptr) {
  if (copy) {
    copy->DeepCopy(orig);
    return nullptr;
  } else if (orig) {
    auto newCopy = vtkSmartPointer<vtkAbstractArray>::Take(orig->NewInstance());
    newCopy->DeepCopy(orig);
    return newCopy;
  } else {
    return nullptr;
  }
}

vtkSmartPointer<vtkAbstractArray> getArrCopy(
    vtkFieldData *field, const std::string &arrayName, int *arrayIdx,
    vtkAbstractArray *array = nullptr) {
  auto orig = arrayIdx ? field->GetAbstractArray(arrayName.c_str(), *arrayIdx)
                       : field->GetAbstractArray(arrayName.c_str());
  return copyOrCreateCopy(orig, array);
}

vtkSmartPointer<vtkAbstractArray> getArrCopy(
    vtkFieldData *field, int arrayIdx, vtkAbstractArray *array = nullptr) {
  auto orig = field->GetAbstractArray(arrayIdx);
  return copyOrCreateCopy(orig, array);
}

}  // namespace

void geoMeshBase::getPointDataArrayCopy(const std::string &arrayName,
                                        vtkAbstractArray *array,
                                        int *arrayIdx) const {
  getArrCopy(_geoMesh.mesh->GetPointData(), arrayName, arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getPointDataArrayCopy(
    const std::string &arrayName, int *arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetPointData(), arrayName, arrayIdx);
}

void geoMeshBase::getPointDataArrayCopy(int arrayIdx,
                                        vtkAbstractArray *array) const {
  getArrCopy(_geoMesh.mesh->GetPointData(), arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getPointDataArrayCopy(
    int arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetPointData(), arrayIdx);
}

void geoMeshBase::getCellDataArrayCopy(const std::string &arrayName,
                                       vtkAbstractArray *array,
                                       int *arrayIdx) const {
  getArrCopy(_geoMesh.mesh->GetCellData(), arrayName, arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getCellDataArrayCopy(
    const std::string &arrayName, int *arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetCellData(), arrayName, arrayIdx);
}

void geoMeshBase::getCellDataArrayCopy(int arrayIdx,
                                       vtkAbstractArray *array) const {
  getArrCopy(_geoMesh.mesh->GetCellData(), arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getCellDataArrayCopy(
    int arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetCellData(), arrayIdx);
}

void geoMeshBase::getFieldDataArrayCopy(const std::string &arrayName,
                                        vtkAbstractArray *array,
                                        int *arrayIdx) const {
  getArrCopy(_geoMesh.mesh->GetFieldData(), arrayName, arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getFieldDataArrayCopy(
    const std::string &arrayName, int *arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetFieldData(), arrayName, arrayIdx);
}

void geoMeshBase::getFieldDataArrayCopy(int arrayIdx,
                                        vtkAbstractArray *array) const {
  getArrCopy(_geoMesh.mesh->GetFieldData(), arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getFieldDataArrayCopy(
    int arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetFieldData(), arrayIdx);
}

// void geoMeshBase::update() {
//  /*
//  if (_stale.GEO_MESH_LINK && !_stale.GEOMETRY && !_stale.MESH) {
//    // Need to loop through all entities in the geometry and all points and
//    // cells in the mesh and set the entity ID.
//    std::vector<GEntity *> entities;
//    _GModel.getEntities(entities, 1);
//    // entities[0]->containsPoint()
//
//    for (vtkIdType ptId = 0; ptId != _mesh->GetNumberOfPoints(); ++ptId) {
//      std::array<double, 3> pt{};
//      _mesh->GetPoint(ptId, pt.data());
//    }
//  }
//  */
//}

geoMeshBase::SideSet::SideSet(vtkPolyData *sideSet, vtkIntArray *geoEnt,
                              vtkIdTypeArray *origCell, vtkIntArray *cellFace,
                              vtkStringArray *setNames)
    : sides(sideSet) {
  assert(geoEnt || sideSet->GetCellData()->HasArray(GEO_ENT_NAME));
  if (geoEnt) { setGeoEntArr(geoEnt); }
  if (origCell) { setOrigCellArr(origCell); }
  if (cellFace) { setCellFaceArr(cellFace); }
  if (setNames) { setSideSetNamesArr(setNames); }
}

geoMeshBase::SideSet::SideSet(vtkPolyData *sides) : SideSet(sides, nullptr) {}

vtkSmartPointer<vtkIntArray> geoMeshBase::SideSet::getGeoEntArr() const {
  return sides ? vtkIntArray::FastDownCast(
                     sides->GetCellData()->GetArray(GEO_ENT_NAME))
               : nullptr;
}

void geoMeshBase::SideSet::setGeoEntArr(vtkIntArray *arr) {
  assert(arr);
  arr->SetName(GEO_ENT_NAME);
  sides->GetCellData()->AddArray(arr);
}

vtkSmartPointer<vtkIdTypeArray> geoMeshBase::SideSet::getOrigCellArr() const {
  return sides ? vtkIdTypeArray::FastDownCast(
                     sides->GetCellData()->GetAbstractArray(ORIG_CELL_NAME))
               : nullptr;
}

void geoMeshBase::SideSet::setOrigCellArr(vtkIdTypeArray *arr) {
  if (arr) {
    assert(arr->GetNumberOfTuples() == sides->GetNumberOfCells() &&
           arr->GetNumberOfComponents() == 2);
    arr->SetName(ORIG_CELL_NAME);
    sides->GetCellData()->AddArray(arr);
  } else {
    sides->GetCellData()->RemoveArray(ORIG_CELL_NAME);
  }
}

vtkSmartPointer<vtkIntArray> geoMeshBase::SideSet::getCellFaceArr() const {
  return sides ? vtkIntArray::FastDownCast(
                     sides->GetCellData()->GetAbstractArray(CELL_FACE_NAME))
               : nullptr;
}

void geoMeshBase::SideSet::setCellFaceArr(vtkIntArray *arr) {
  if (arr) {
    assert(arr->GetNumberOfTuples() == sides->GetNumberOfCells() &&
           arr->GetNumberOfComponents() == 2);
    arr->SetName(CELL_FACE_NAME);
    sides->GetCellData()->AddArray(arr);
  } else {
    sides->GetCellData()->RemoveArray(CELL_FACE_NAME);
  }
}

vtkSmartPointer<vtkStringArray> geoMeshBase::SideSet::getSideSetNames() const {
  return sides ? vtkStringArray::SafeDownCast(
                     sides->GetFieldData()->GetAbstractArray(NAME_ARR_NAME))
               : nullptr;
}

void geoMeshBase::SideSet::setSideSetNamesArr(vtkStringArray *arr) {
  if (arr) {
    assert(arr->GetNumberOfComponents() == 1);
    arr->SetName(NAME_ARR_NAME);
    sides->GetFieldData()->AddArray(arr);
  } else {
    sides->GetFieldData()->RemoveArray(NAME_ARR_NAME);
  }
}

void geoMeshBase::GeoMesh::findSide2OrigCell() {
  if (sideSet.sides && !sideSet.getOrigCellArr() && !sideSet.getCellFaceArr()) {
    vtkIdType numSides = sideSet.sides->GetNumberOfCells();
    auto entArr = sideSet.getGeoEntArr();
    vtkNew<vtkIdTypeArray> origCellArr;
    origCellArr->SetNumberOfComponents(2);
    origCellArr->SetNumberOfTuples(numSides);
    origCellArr->FillTypedComponent(0, -1);
    origCellArr->FillTypedComponent(1, -1);
    vtkNew<vtkIntArray> cellFaceArr;
    cellFaceArr->SetNumberOfComponents(2);
    cellFaceArr->SetNumberOfTuples(numSides);
    cellFaceArr->FillTypedComponent(0, -1);
    cellFaceArr->FillTypedComponent(1, -1);
    // Using points->cell, intersect the set of cells that have the same points
    // as the side cell
    mesh->BuildLinks();
    auto pt2Cell = mesh->GetCellLinks();
    for (vtkIdType i = 0; i < sideSet.sides->GetNumberOfCells(); ++i) {
      auto sidePoints = sideSet.sides->GetCell(i)->GetPointIds();
      auto cellList = pt2Cell->GetLink(sidePoints->GetId(0));
      std::vector<vtkIdType> candidates{cellList.cells,
                                        cellList.cells + cellList.ncells};
      std::sort(candidates.begin(), candidates.end());
      for (auto j = 1; j < sidePoints->GetNumberOfIds(); ++j) {
        std::vector<vtkIdType> newCandidates{};
        cellList = pt2Cell->GetLink(sidePoints->GetId(j));
        for (auto iterCell = cellList.cells;
             iterCell < cellList.cells + cellList.ncells; ++iterCell) {
          auto findIter =
              std::lower_bound(candidates.begin(), candidates.end(), *iterCell);
          if (findIter != candidates.end() && *iterCell == *findIter) {
            newCandidates.emplace_back(*iterCell);
          }
        }
        candidates = std::move(newCandidates);
        std::sort(candidates.begin(), candidates.end());
      }
      // For each candidate cell, check each edge/face for a match
      bool foundMatch = false;
      if (candidates.size() == 1) {
        // If there is only one candidate (side is on boundary, not material
        // interface), check by sorting points - note no check for orientation
        std::vector<vtkIdType> sortedSidePoints(sidePoints->GetNumberOfIds());
        auto pointsBegin = sidePoints->GetPointer(0);
        std::partial_sort_copy(
            pointsBegin, pointsBegin + sidePoints->GetNumberOfIds(),
            sortedSidePoints.begin(), sortedSidePoints.end());
        auto candOrigCell = mesh->GetCell(candidates.at(0));
        auto dim = candOrigCell->GetCellDimension();
        for (int j = 0; j < (dim == 3 ? candOrigCell->GetNumberOfFaces()
                                      : candOrigCell->GetNumberOfEdges());
             ++j) {
          auto candSide =
              dim == 3 ? candOrigCell->GetFace(j) : candOrigCell->GetEdge(j);
          std::vector<vtkIdType> sortedCandPoints(sortedSidePoints.size());
          auto candPointsBegin = candSide->GetPointIds()->GetPointer(0);
          // Note use of sidePoints->GetNumberOfIds() in case, eg, candSide is
          // a vtkQuadraticTriangle, but side is a vtkTriangle
          auto candPointsEnd = candPointsBegin + sidePoints->GetNumberOfIds();
          std::partial_sort_copy(candPointsBegin, candPointsEnd,
                                 sortedCandPoints.begin(),
                                 sortedCandPoints.end());
          if (sortedSidePoints == sortedCandPoints) {
            foundMatch = true;
            origCellArr->SetTypedComponent(i, 0, candidates.at(0));
            cellFaceArr->SetTypedComponent(i, 0, j);
            break;
          }
        }
      } else {
        // If there are two candidates, match the orientation of the side
        assert(candidates.size() == 2);
        for (const auto &candIdx : candidates) {
          auto candOrigCell = mesh->GetCell(candIdx);
          auto dim = candOrigCell->GetCellDimension();
          bool foundSide = false;
          for (int j = 0; j < (dim == 3 ? candOrigCell->GetNumberOfFaces()
                                        : candOrigCell->GetNumberOfEdges()) &&
                          !foundSide;
               ++j) {
            auto candSide =
                dim == 3 ? candOrigCell->GetFace(j) : candOrigCell->GetEdge(j);
            auto candPointsBegin = candSide->GetPointIds()->GetPointer(0);
            // Note use of sidePoints->GetNumberOfIds() in case, eg, candSide is
            // a vtkQuadraticTriangle, but side is a vtkTriangle
            auto candPointsEnd = candPointsBegin + sidePoints->GetNumberOfIds();
            auto findFirstPoint =
                std::find(candPointsBegin, candPointsEnd, sidePoints->GetId(0));
            if (findFirstPoint != candPointsEnd) {
              std::vector<vtkIdType> rotatedCandPoints(
                  sidePoints->GetNumberOfIds());
              std::rotate_copy(candPointsBegin, findFirstPoint, candPointsEnd,
                               rotatedCandPoints.begin());
              if (std::equal(rotatedCandPoints.begin(), rotatedCandPoints.end(),
                             sidePoints->GetPointer(0))) {
                foundMatch = true;
                foundSide = true;
                origCellArr->SetTypedComponent(i, 0, candIdx);
                cellFaceArr->SetTypedComponent(i, 0, j);
              } else if (std::equal(rotatedCandPoints.rbegin(),
                                    rotatedCandPoints.rend(),
                                    sidePoints->GetPointer(0))) {
                foundSide = true;
                origCellArr->SetTypedComponent(i, 1, candIdx);
                cellFaceArr->SetTypedComponent(i, 1, j);
              }
            }
          }
          if (foundMatch) { break; }
        }
      }
      if (!foundMatch) {
        std::cerr << "No cell in mesh found to correspond to side set cell "
                  << i << ".\n";
      }
    }
    sideSet.setOrigCellArr(origCellArr);
    sideSet.setCellFaceArr(cellFaceArr);
  }
}

int geoMeshBase::GeoMesh::getDimension() const {
  if (this->mesh->GetNumberOfCells() > 0) {
    return this->mesh->GetCell(0)->GetCellDimension();
  } else {
    return -1;
  }
}

}  // namespace MSH
}  // namespace NEM
