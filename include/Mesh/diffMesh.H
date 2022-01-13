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
#ifndef NEMOSYS_DIFFMESH_H_
#define NEMOSYS_DIFFMESH_H_

#include "nemosys_export.h"
#include "Mesh/geoMeshBase.H"

namespace NEM {
namespace MSH {

/**
 * Compare two @c geoMeshBase objects.
 * @param gmb1 Reference @c geoMeshBase
 * @param gmb2 Other @c geoMeshBase
 * @param floor If a floating point value is less than or equal to @p floor, it
 * is treated as 0. Should be non-negative
 * @param relTol Relative tolerance for comparing floating point values
 * @param numCellsTol Tolerance for number of cells, if positive. If non-positive,
 * exact match required for number of cells and cell data compared at each cell.
 * @param numPointsTol Tolerance for number of points, if positive. If
 * non-positive, exact match required for number of points and point data
 * compared at each cell.
 * @return 0 if @p gmb1 and @p gmb2 are similar; 1 otherwise
 */
NEMOSYS_EXPORT int diffMesh(geoMeshBase *gmb1, geoMeshBase *gmb2,
                            double floor = 1e-9, double relTol = 1e-6,
                            double numCellsTol = -1.,
                            double numPointsTol = -1.);

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_DIFFMESH_H_
