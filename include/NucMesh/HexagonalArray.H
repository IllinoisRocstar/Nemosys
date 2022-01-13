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
#ifndef NEMOSYS_HEXAGONALARRAY_H_
#define NEMOSYS_HEXAGONALARRAY_H_

#include "nemosys_export.h"
#include "NucMesh/ShapesArray.H"

#include <array>

namespace NEM {
namespace NUCMESH {

class NEMOSYS_EXPORT HexagonalArray : public ShapesArray {
 public:
  HexagonalArray(std::size_t numRadii, double deltaRadius,
                 const std::array<double, 3> &center = {0, 0, 0});

  const std::size_t &getPatternRowCol(int row, int col) const;
  const std::size_t &getPatternCoordCenter(int right, int rightUp) const;

  void setPatternRowCol(int row, int col, std::size_t patternKey);
  /**
   * The shape will be translated by <tt> right * [getGridDistance(), 0, 0] +
   * rightUp * [getGridDistance() * cos(60), getGridDistance() * sin(60), 0] +
   * getCenter() </tt>.
   */
  void setPatternCoordCenter(int right, int rightUp, std::size_t patternKey);

  std::size_t getNumRadii() const { return numRadii_; }
  double getGridDistance() const { return delta_; }
  void setGridDistance(double gridDistance) { delta_ = gridDistance; }

 protected:
  NEM::GEO::GeoManager createGeo() const override;

 private:
  double delta_;
  std::size_t numRadii_;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_HEXAGONALARRAY_H_
