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
#include "NucMesh/RectangularArray.H"

#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

namespace NEM {
namespace NUCMESH {

RectangularArray::RectangularArray(const std::array<std::size_t, 2> &gridDims,
                                   const std::array<double, 2> &deltaGrid,
                                   const std::array<double, 3> &center)
    : ShapesArray(center, gridDims[0] * gridDims[1]),
      dims_(gridDims),
      delta_(deltaGrid) {}

const std::size_t &RectangularArray::getPattern(std::size_t x,
                                                std::size_t y) const {
  return this->ShapesArray::getPattern(x * dims_[1] + y);
}

void RectangularArray::setPattern(std::size_t x, std::size_t y,
                                  std::size_t patternKey) {
  this->ShapesArray::setPattern(x * dims_[1] + y, patternKey);
}

NEM::GEO::GeoManager RectangularArray::createGeo() const {
  const auto maxY = static_cast<int>(this->dims_[1] - 1);
  std::array<int, 2> coord{-static_cast<int>(this->dims_[0] - 1), -maxY};
  return this->createGeoImpl([this, maxY,
                              &coord](NEM::GEO::GeoManager *const inp) {
    if (inp) {
      auto &center = this->getCenter();
      gp_Trsf transformation{};
      transformation.SetTranslation(
          gp_Vec{center[0] + coord[0] * this->delta_[0] / 2.,
                 center[1] + coord[1] * this->delta_[1] / 2., center[2]});
      *inp = ShapesArray::basicTransformation(transformation, std::move(*inp));
    }
    coord[1] += 2;
    if (coord[1] > maxY) {
      coord[1] = -maxY;
      coord[0] += 2;
    }
  });
}

}  // namespace NUCMESH
}  // namespace NEM