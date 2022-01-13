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
#ifndef NEMOSYS_PACKMESHDRIVER_H_
#define NEMOSYS_PACKMESHDRIVER_H_

#include "Drivers/NemDriver.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace NEM {
namespace DRV {

/**
  @brief This class drives the flow of pack meshing by using parameters provided
  by users. It supports two different meshing pipelines. First being hexahedral
  conformal meshing and second is tetrahedral periodic as well as non-priodic
  meshing.

  Tetrahedral periodic meshing pipeline operates within rocPack class and uses
  Netgen mesh optimization for final mesh quality improvement if needed.
  ## Pack Meshing Workflow for Tetrahedral Mesh
  > * Geometry generation using RocPack output file.
  > * Enforcing periodicity in geometry for periodic mesh
  > * Mapping periodic surfaces and points with their counterparts.
  > * Meshing the geometry tetrahedrally while enforcing periodic boundary
  boundary conditions.
  > * Exporting final mesh with or without optimization.

**/
class NEMOSYS_EXPORT PackMeshDriver : public NemDriver {
  JSONCONS_TYPE_TRAITS_FRIEND
 protected:
  PackMeshDriver() = default;

 private:
  static constexpr const char *programType = "Pack Mesh Generation";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_PACKMESHDRIVER_H_
