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
#ifndef NEMOSYS_SHAPEDATA_H_
#define NEMOSYS_SHAPEDATA_H_

#include "nemosys_export.h"

#include <memory>
#include <vector>

class TopoDS_Shape;
template <typename T>
class NCollection_List;
using TopTools_ListOfShape = NCollection_List<TopoDS_Shape>;

namespace NEM {
namespace GEO {

class GeoManager;

/**
 * Data stored in an instance of @c GeoManager
 */
class NEMOSYS_EXPORT ShapeData {
 public:
  virtual ~ShapeData() = default;

  /**
   * @brief Method that is called when a geometry operation generates new shapes
   * @details Default implementation does nothing
   * @param[in] shape Argument to operation; shape existing in @p geoMetadata
   * @param[in] generatedShapes Shapes provided by a call to @c
   * op.Generated(shape) for some geometry operation
   * @param[in,out] geoMetadata Geometry-to-data map to alter
   */
  virtual void updateGenerated(const TopoDS_Shape &shape,
                               const TopTools_ListOfShape &generatedShapes,
                               GeoManager &geoMetadata) const {};

  /**
   * @brief Method that is called when a geometry operation modifies @p shape
   * @param[in] shape Argument to operation; shape existing in @p geoMetadata
   * @param[in] modifiedShapes Shapes provided by a call to @c
   * op.Modified(shape) for some geometry operation
   * @param[in,out] geoMetadata Geometry-to-data map to alter
   */
  virtual void updateModified(const TopoDS_Shape &shape,
                              const TopTools_ListOfShape &modifiedShapes,
                              GeoManager &geoMetadata) const = 0;

  /**
   * @brief Method that is called when @p geoMetadata removes @p shape
   * @details Default implementation does nothing
   * @param[in] shape Shape to be removed
   * @param[in,out] geoMetadata Geometry-to-data map to alter
   */
  virtual void updateDeleted(const TopoDS_Shape &shape,
                             GeoManager &geoMetadata) {};

 protected:
  ShapeData() = default;
  ShapeData(const ShapeData &) = default;
  ShapeData(ShapeData &&) = default;
  ShapeData &operator=(const ShapeData &) = default;
  ShapeData &operator=(ShapeData &&) = default;
};

}  // namespace GEO
}  // namespace NEM

#endif  // NEMOSYS_SHAPEDATA_H_
