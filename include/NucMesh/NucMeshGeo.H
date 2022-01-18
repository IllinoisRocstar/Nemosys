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
#ifndef NEMOSYS_NUCMESHGEO_H_
#define NEMOSYS_NUCMESHGEO_H_

#include "nemosys_export.h"
#include "Geometry/GeoManager.H"

#include <memory>
#include <utility>

class SMESH_Mesh;
class SMESH_Gen;

namespace NEM {
namespace NUCMESH {

class NEMOSYS_EXPORT NucMeshGeo : public NEM::GEO::GeoManager {
  using NEM::GEO::GeoManager::GeoManager;
 public:
  /**
   * @brief Construct a mesh using the shapes in the map and their @c ShapeData.
   * @details Note there is a default edge hypothesis
   * (StdMeshers_AutomaticLength), edge algorithm (StdMeshers_Regular_1D), and
   * face algorithm (StdMeshers_MEFISTO_2D). Geometry and mesh regenerated on
   * each call.
   * @return Generated mesh
   */
  std::unique_ptr<SMESH_Mesh> computeMesh(SMESH_Gen &generator);
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_NUCMESHGEO_H_
