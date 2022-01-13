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
#ifndef NEMOSYS_RECTANGULARARRAY_H_
#define NEMOSYS_RECTANGULARARRAY_H_

#include "nemosys_export.h"
#include "NucMesh/ShapesArray.H"

#include <array>

namespace NEM {
namespace NUCMESH {

class NEMOSYS_EXPORT RectangularArray : public ShapesArray {
 public:
  RectangularArray(const std::array<std::size_t, 2> &gridDims,
                   const std::array<double, 2> &deltaGrid,
                   const std::array<double, 3> &center = {0, 0, 0});

  const std::size_t &getPattern(std::size_t x, std::size_t y) const;
  /**
   * The shape at <tt> [x, y] </tt> will be translated by <tt> [(x -
   * (getGridDims()[0] - 1) / 2) * getGridDistance()[0], (y - (getGridDims()[1]
   * - 1) / 2) * getGridDistance()[1], 0] + getCenter() </tt>.
   */
  void setPattern(std::size_t x, std::size_t y, std::size_t patternKey);

  const std::array<std::size_t, 2> &getGridDims() const { return dims_; }
  const std::array<double, 2> &getGridDistance() const { return delta_; }
  void setGridDistance(const std::array<double, 2> &gridDistance) {
    delta_ = gridDistance;
  }

 protected:
  NEM::GEO::GeoManager createGeo() const override;

 private:
  std::array<std::size_t, 2> dims_;
  std::array<double, 2> delta_;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_RECTANGULARARRAY_H_
