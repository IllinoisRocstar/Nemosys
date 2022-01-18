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
#ifndef NEMOSYS_NUCMESHDRIVER_H_
#define NEMOSYS_NUCMESHDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <vtkSmartPointer.h>

#include "Services/NucMeshSrv.H"

namespace NEM {
namespace DRV {

/**
 * @brief Class to create 2d geometry from polygons and circles and mesh it.
 */
class NEMOSYS_EXPORT NucMeshDriver : public NemDriver {
 public:
  using Files = DriverOutFile;

  using Opts = NEM::SRV::NucMeshConf;

  /**
   * NucMeshDriver ctor
   * @param file Output file. Should end in ".vtu"
   * @param opts Specify shapes and extrusion.
   */
  explicit NucMeshDriver(Files file, Opts opts);

  /**
   * Get the output file name
   * @return output file object containing file name
   */
  const Files &getFiles() const;
  /**
   * Set the output file for the driver
   * @param files NEM::DRV::NucMeshDriver::Files object containing output file
   */
  void setFiles(Files files);

  /**
   * Get the options set on a driver object
   * @return geometry, meshing options, and extrusion steps
   */
  const Opts &getOpts() const;
  /**
   * Set the geometry/mesh and extrusion to execute
   * @param opts geometry/mesh options
   */
  void setOpts(Opts opts);

  /**
   * Construct the mesh
   * @return geoMeshBase with 1d and 2d cells (and 3d if any extrusion)
   */
  vtkSmartPointer<NEM::MSH::geoMeshBase> draw() const;

  /**
   * Construct the mesh and write to file.
   */
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  NucMeshDriver();
  Files file_;
  Opts opts_;

  static constexpr const char *programType = "NucMesh Generation";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_NUCMESHDRIVER_H_
