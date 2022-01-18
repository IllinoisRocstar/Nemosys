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
#include "NucMesh/ShapesArray.H"

#include <memory>

#include <BRepBuilderAPI_Transform.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Compound.hxx>

namespace NEM {
namespace NUCMESH {

ShapesArray::ShapesArray(const std::array<double, 3> &center,
                         std::size_t numPatternShapes)
    : ShapeBase(center), patternShapes_(), pattern_(numPatternShapes) {}

std::size_t ShapesArray::getNumPatternShapes() const {
  return patternShapes_.size();
}

const ShapeBase *ShapesArray::getPatternShape(std::size_t idx) const {
  if (idx < patternShapes_.size()) {
    return patternShapes_.at(idx).get();
  } else {
    return nullptr;
  }
}

void ShapesArray::setPatternShape(std::size_t idx,
                                  const std::shared_ptr<ShapeBase> &shape) {
  if (idx >= patternShapes_.size()) {
    patternShapes_.resize(idx + 1);
  }
  patternShapes_.at(idx) = shape;
}

void ShapesArray::fillPattern(std::size_t idx) {
  for (auto &pattern : pattern_) {
    pattern = idx;
  }
}

std::size_t ShapesArray::getPatternSize() const { return pattern_.size(); }

const std::size_t &ShapesArray::getPattern(std::size_t idx) const {
  return pattern_.at(idx);
}

void ShapesArray::setPattern(std::size_t idx, std::size_t patternKey) {
  pattern_.at(idx) = patternKey;
}

NEM::GEO::GeoManager ShapesArray::basicTransformation(
    const gp_Trsf &transformation, NEM::GEO::GeoManager &&geoMetadata) {
  auto compound = geoMetadata.buildCompound();
  BRepBuilderAPI_Transform transformer{transformation};
  transformer.Perform(compound);
  NEM::GEO::GeoManager output(geoMetadata.getDim());
  static constexpr std::array<TopAbs_ShapeEnum, 4> shapeTypes{
      TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX};
  for (auto &shapeType : shapeTypes) {
    for (TopExp_Explorer explorer{compound, shapeType}; explorer.More();
         explorer.Next()) {
      auto &oldSubshape = explorer.Current();
      if (auto old_metadata = geoMetadata.get(oldSubshape)) {
        output.insert(transformer.ModifiedShape(oldSubshape),
                      std::move(*old_metadata));
      }
    }
  }
  return output;
}

}  // namespace NUCMESH
}  // namespace NEM