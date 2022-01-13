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
#include "MeshOperation/dataSetRegionBoundaryFilter.H"

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>

namespace NEM {
namespace MSH {

namespace {

struct BoundaryCell {
  BoundaryCell(vtkIdType source, int side, int region)
      : sources{source, -1}, sides{side, -1}, regions{region, -1} {}

  void setTwin(vtkIdType otherSource, int otherSide, int otherRegion) {
    this->sources[1] = otherSource;
    this->sides[1] = otherSide;
    this->regions[1] = otherRegion;
    if (otherRegion < this->regions[0]) {
      std::swap(this->sources[0], this->sources[1]);
      std::swap(this->sides[0], this->sides[1]);
      std::swap(this->regions[0], this->regions[1]);
    }
  }

  std::array<vtkIdType, 2> sources;
  std::array<int, 2> sides;
  std::array<int, 2> regions;
};

class BoundaryCellSet {
 public:
  using data_type = std::vector<
      std::vector<std::pair<std::array<vtkIdType, 2>, BoundaryCell>>>;
  explicit BoundaryCellSet(data_type::size_type size)
      : data(size), total_size(0) {}

  class const_iterator
      : public std::iterator<std::forward_iterator_tag, const BoundaryCell> {
   private:
    const_iterator(data_type::const_iterator outerIter,
                   data_type::const_iterator endOuterIter,
                   data_type::value_type::const_iterator innerIter)
        : outerIter(outerIter),
          endOuterIter(endOuterIter),
          innerIter(innerIter){};

   public:
    reference operator*() { return innerIter->second; }
    pointer operator->() { return &innerIter->second; }
    const_iterator &operator++() {
      ++innerIter;
      if (innerIter == outerIter->end()) {
        do {
          ++outerIter;
          if (outerIter == endOuterIter) {
            return *this;
          }
        } while (outerIter->empty());
        innerIter = outerIter->begin();
      }
      return *this;
    }
    const_iterator operator++(int) {
      auto temp = *this;
      ++(*this);
      return temp;
    }
    friend bool operator==(const const_iterator &l, const const_iterator &r) {
      return (l.outerIter == r.outerIter && l.innerIter == r.innerIter) ||
             (l.outerIter == l.endOuterIter && r.outerIter == r.endOuterIter);
    }
    friend bool operator!=(const const_iterator &l, const const_iterator &r) {
      return !(l == r);
    }
    friend BoundaryCellSet;

   private:
    data_type::const_iterator outerIter;
    const data_type::const_iterator endOuterIter;
    data_type::value_type::const_iterator innerIter;
  };

  const_iterator begin() const {
    auto beginOuter = data.begin();
    const auto endOuter = data.end();
    while (beginOuter->empty()) {
      ++beginOuter;
      if (beginOuter == endOuter) {
        return {beginOuter, endOuter, {}};
      }
    }
    return {beginOuter, endOuter, beginOuter->begin()};
  }

  const_iterator end() const {
    auto endIter = data.end();
    return {endIter, endIter, {}};
  }

  void insert(vtkIdType source, int side, vtkIdList *points, int numPoints,
              int region) {
    auto pointsBegin = points->GetPointer(0);
    auto numPointsKey = std::min(numPoints, 3);
    auto pointsEnd = pointsBegin + numPoints;
    auto pointsMiddle = pointsBegin + numPointsKey;
    // We only need the least 2/3 point ids.
    std::partial_sort(pointsBegin, pointsMiddle, pointsEnd);
    insertImpl(data[pointsBegin[1]],
               {pointsBegin[0], numPointsKey < 3 ? -1 : pointsBegin[2]}, source,
               side, region);
  }

  void insertTri(vtkIdType source, int face, vtkIdType point0, vtkIdType point1,
                 vtkIdType point2, int region) {
    if (point0 > point1) std::swap(point0, point1);
    if (point0 > point2) std::swap(point0, point2);
    if (point1 > point2) std::swap(point1, point2);
    insertImpl(this->data[point1], {point0, point2}, source, face, region);
  }

  void insertQuad(vtkIdType source, int face, vtkIdType point0,
                  vtkIdType point1, vtkIdType point2, vtkIdType point3,
                  int region) {
    if (point0 > point1) std::swap(point0, point1);
    if (point2 > point3) std::swap(point2, point3);
    if (point0 > point2) std::swap(point0, point2);
    if (point1 > point3) std::swap(point1, point3);
    if (point1 > point2) std::swap(point1, point2);
    insertImpl(this->data[point1], {point0, point2}, source, face, region);
  }

  // These methods are faster than repeated calls to cell->GetFace() then
  // this->insert.
  void insertTetFaces(vtkIdType source, const vtkIdType *points, int region) {
    insertTri(source, 0, points[0], points[1], points[3], region);
    insertTri(source, 1, points[1], points[2], points[3], region);
    insertTri(source, 2, points[2], points[0], points[3], region);
    insertTri(source, 3, points[0], points[2], points[1], region);
  }

  void insertWedgeFaces(vtkIdType source, const vtkIdType *points, int region) {
    insertTri(source, 0, points[0], points[1], points[2], region);
    insertTri(source, 1, points[3], points[5], points[4], region);
    insertQuad(source, 2, points[0], points[3], points[4], points[1], region);
    insertQuad(source, 3, points[1], points[4], points[5], points[2], region);
    insertQuad(source, 4, points[2], points[5], points[3], points[0], region);
  }

  void insertHexFaces(vtkIdType source, const vtkIdType *points, int region) {
    insertQuad(source, 0, points[0], points[4], points[7], points[3], region);
    insertQuad(source, 1, points[1], points[2], points[6], points[5], region);
    insertQuad(source, 2, points[0], points[1], points[5], points[4], region);
    insertQuad(source, 3, points[3], points[7], points[6], points[2], region);
    insertQuad(source, 4, points[0], points[3], points[2], points[1], region);
    insertQuad(source, 5, points[4], points[5], points[6], points[7], region);
  }

  data_type::value_type::size_type size() const { return total_size; }

 private:
  /**
   * Edge/face with sorted point ids (a, b, c, ...) is located at some index i
   * in data[b], with data[b][i].first == [a, c] (for edges, third point id
   * treated as -1). Using the second-least point as hash results in fewer
   * collisions (using sum is even better as a hash, but slower to compute).
   */
  data_type data;
  data_type::value_type::size_type total_size;

  void insertImpl(std::decay<decltype(data)>::type::value_type &bucket,
                  std::array<vtkIdType, 2> key, int source, int side,
                  int region) {
    auto matchFind = std::find_if(
        bucket.begin(), bucket.end(),
        [&key](const std::decay<decltype(bucket)>::type::value_type &x) {
          return x.first == key;
        });
    if (matchFind != bucket.end()) {
      if (matchFind->second.regions[0] == region) {
        // Remove from bucket
        *matchFind = bucket.back();
        bucket.pop_back();
        total_size -= 1;
      } else {
        // Set twin information
        matchFind->second.setTwin(source, side, region);
      }
    } else {
      // Add to bucket
      bucket.emplace_back(std::piecewise_construct, std::forward_as_tuple(key),
                          std::forward_as_tuple(source, side, region));
      total_size += 1;
    }
  }
};

}  // namespace

vtkStandardNewMacro(dataSetRegionBoundaryFilter)

int dataSetRegionBoundaryFilter::FillInputPortInformation(
    int port, vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

int dataSetRegionBoundaryFilter::RequestData(
    vtkInformation *request, vtkInformationVector **inputVector,
    vtkInformationVector *outputVector) {
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  auto input =
      vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  auto output =
      vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (input->CheckAttributes()) {
    return 0;
  }

  const vtkIdType numCells = input->GetNumberOfCells();
  if (numCells == 0) {
    vtkDebugMacro(<< "Number of cells is zero, no data to process.");
    return 1;
  }
  vtkNew<vtkCellArray> newCells{};
  // Approximate
  newCells->Allocate(numCells);

  BoundaryCellSet boundaryCellSet(input->GetNumberOfPoints());
  auto regionArr = vtkArrayDownCast<vtkIntArray>(
      input->GetCellData()->GetAbstractArray(this->MaterialArrayName.c_str()));
  for (vtkIdType i = 0; i < numCells; ++i) {
    auto cell = input->GetCell(i);
    auto region = regionArr->GetTypedComponent(i, 0);

    if (cell->GetCellDimension() != this->GetDimension()) {
      continue;
    }

    if (this->GetDimension() == 2) {
      for (int j = 0; j < cell->GetNumberOfEdges(); ++j) {
        auto edge = cell->GetEdge(j);
        boundaryCellSet.insert(i, j, edge->GetPointIds(), 2, region);
      }
    } else {  // this->GetDimension() == 3
      switch (static_cast<VTKCellType>(cell->GetCellType())) {
        case VTK_TETRA:
        case VTK_QUADRATIC_TETRA:
        case VTK_HIGHER_ORDER_TETRAHEDRON:
        case VTK_LAGRANGE_TETRAHEDRON: {
          auto points = cell->GetPointIds()->GetPointer(0);
          boundaryCellSet.insertTetFaces(i, points, region);
          break;
        }
        case VTK_HEXAHEDRON:
        case VTK_QUADRATIC_HEXAHEDRON:
        case VTK_TRIQUADRATIC_HEXAHEDRON:
        case VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON:
        case VTK_HIGHER_ORDER_HEXAHEDRON:
        case VTK_LAGRANGE_HEXAHEDRON: {
          auto points = cell->GetPointIds()->GetPointer(0);
          boundaryCellSet.insertHexFaces(i, points, region);
          break;
        }
        case VTK_WEDGE:
        case VTK_QUADRATIC_WEDGE:
        case VTK_QUADRATIC_LINEAR_WEDGE:
        case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
        case VTK_HIGHER_ORDER_WEDGE:
        case VTK_LAGRANGE_WEDGE: {
          auto points = cell->GetPointIds()->GetPointer(0);
          boundaryCellSet.insertWedgeFaces(i, points, region);
          break;
        }
        default: {
          for (int j = 0; j < cell->GetNumberOfFaces(); ++j) {
            auto face = cell->GetFace(j);
            // Treat non-linear faces as linear ones by only looking at
            // beginning of point id array. Number of edges should be equal to
            // number of vertices even for non-linear faces.
            boundaryCellSet.insert(i, j, face->GetPointIds(),
                                   face->GetNumberOfEdges(), region);
          }
          break;
        }
      }
    }
  }

  vtkSmartPointer<vtkIdTypeArray> origCellId{};
  vtkSmartPointer<vtkIntArray> cellSideId{}, surfRegionId{};
  if (!this->GetOrigCellArrayName().empty()) {
    origCellId = vtkSmartPointer<vtkIdTypeArray>::New();
    origCellId->SetName(this->GetOrigCellArrayName().c_str());
    origCellId->SetNumberOfComponents(2);
    origCellId->SetNumberOfTuples(boundaryCellSet.size());
  }
  if (!this->GetCellSideArrayName().empty()) {
    cellSideId = vtkSmartPointer<vtkIntArray>::New();
    cellSideId->SetName(this->GetCellSideArrayName().c_str());
    cellSideId->SetNumberOfComponents(2);
    cellSideId->SetNumberOfTuples(boundaryCellSet.size());
  }
  if (!this->GetBoundaryRegionArrayName().empty()) {
    surfRegionId = vtkSmartPointer<vtkIntArray>::New();
    surfRegionId->SetName(this->GetBoundaryRegionArrayName().c_str());
    surfRegionId->SetNumberOfTuples(boundaryCellSet.size());
  }

  // Map from regionArr to surfRegionId.
  std::map<std::array<int, 2>, int> oldRegions2SurfRegions{};
  for (const auto &boundaryCell : boundaryCellSet) {
    auto origCellToInsert = input->GetCell(boundaryCell.sources[0]);
    auto sideToInsert = this->GetDimension() == 2
                            ? origCellToInsert->GetEdge(boundaryCell.sides[0])
                            : origCellToInsert->GetFace(boundaryCell.sides[0]);
    auto sideIdx = newCells->InsertNextCell(
        this->GetDimension() == 2 ? 2 : sideToInsert->GetNumberOfEdges(),
        sideToInsert->GetPointIds()->GetPointer(0));
    if (origCellId) {
      origCellId->SetTypedTuple(sideIdx, boundaryCell.sources.data());
    }
    if (cellSideId) {
      cellSideId->SetTypedTuple(sideIdx, boundaryCell.sides.data());
    }
    auto surfRegionIter = oldRegions2SurfRegions.emplace(
        boundaryCell.regions, oldRegions2SurfRegions.size());
    if (surfRegionId) {
      surfRegionId->SetTypedComponent(sideIdx, 0, surfRegionIter.first->second);
    }
  }
  newCells->Squeeze();

  if (!this->GetRegionToMaterialArrayName().empty()) {
    vtkNew<vtkIntArray> region2MaterialId;
    region2MaterialId->SetName(this->GetRegionToMaterialArrayName().c_str());
    region2MaterialId->SetNumberOfComponents(2);
    region2MaterialId->SetNumberOfTuples(oldRegions2SurfRegions.size());
    for (const auto &old2newRegion : oldRegions2SurfRegions) {
      auto newRegion = old2newRegion.second;
      region2MaterialId->SetComponent(newRegion, 0, old2newRegion.first[0]);
      region2MaterialId->SetComponent(newRegion, 1, old2newRegion.first[1]);
    }
    output->GetFieldData()->AddArray(region2MaterialId);
  }

  if (auto pointSet = vtkPointSet::SafeDownCast(input)) {
    output->SetPoints(pointSet->GetPoints());
  } else if (auto rectGrid = vtkRectilinearGrid::SafeDownCast(input)) {
    rectGrid->GetPoints(output->GetPoints());
  } else {
    auto outPoints = output->GetPoints();
    outPoints->Initialize();
    outPoints->SetNumberOfPoints(input->GetNumberOfPoints());
    for (vtkIdType i = 0; i < input->GetNumberOfPoints(); ++i) {
      outPoints->SetPoint(i, input->GetPoint(i));
    }
  }
  output->GetPointData()->ShallowCopy(input->GetPointData());

  if (this->GetDimension() == 2) {
    output->SetLines(newCells);
  } else {  // this->GetDimension() == 3
    output->SetPolys(newCells);
  }

  vtkCellData *outputCD = output->GetCellData();
  outputCD->AddArray(origCellId);
  outputCD->AddArray(cellSideId);
  outputCD->AddArray(surfRegionId);

  if (output->CheckAttributes()) {
    return 0;
  } else {
    return 1;
  }
}

}  // namespace MSH
}  // namespace NEM
