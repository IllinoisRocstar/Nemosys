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
#ifndef NEMOSYS_SURFACEPACKMESHDRIVER_H_
#define NEMOSYS_SURFACEPACKMESHDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/PackMesh/PackMeshDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT SurfacePackMeshDriver : public PackMeshDriver {
 public:
  struct NEMOSYS_EXPORT Files {
   public:
    /**
     * @brief
     * @param rocpackFile Input file name (Output from Rocpack)
     * @param outputMeshFile Output file name (.stl, .vtu, .msh are supported)
     */
    Files(std::string rocpackFile, std::string outputMeshFile);
    std::string rocpackFile{};
    std::string outputMeshFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Files() = default;
  };

  enum class NEMOSYS_EXPORT PhysGrpOpts { NONE, MULTI, TWO, PER_SHAPE };

  struct NEMOSYS_EXPORT CustomDomain {
    std::array<double, 3> initial{};
    std::array<double, 3> length{};
  };

  struct NEMOSYS_EXPORT Periodic3DOpts {
    PhysGrpOpts physGrpOptions{PhysGrpOpts::NONE};
    /**
     * @brief Boolean for enabling cohesive elements creation
     */
    bool createCohesive{false};
    /**
     * @brief Boolean for enabling patches.
     */
    bool enablePatches{false};
    /**
     * @brief Boolean for setting periodic geometry
     */
    bool setPeriodicGeo{false};
    /**
     * @brief for users to specify custom domain size.
     */
    jsoncons::optional<CustomDomain> customDomain{};
    /**
     * @brief Parameters for translating mesh from origin.
     */
    std::array<double, 3> transferMesh{};
    /**
     * @brief Element order selected by user
     */
    int elemOrder{1};
  };

  struct NEMOSYS_EXPORT Opts {
   public:
    explicit Opts(jsoncons::optional<Periodic3DOpts> periodic3DOpts);
    /**
     * @brief Mesh size
     */
    jsoncons::optional<double> meshSize{};
    /**
     * @brief Meshing algorithm of choice.
     */
    int meshAlgorithm{1};
    /**
     * @brief Value of scaling between 0 to 1 with 1 being no scaling.
     */
    jsoncons::optional<double> scaleValue{};
    /**
     * @brief Boolean for removal of volumes intersecting boundary
     */
    bool removeBoundaryPacks{false};
    /**
     * @brief Enable default output formats
     */
    bool enableDefaultOut{false};
    /**
     * @brief Boolean for preserving pack sizes instead of packing fraction
     */
    bool preserveSize{false};
    /**
     * @brief Refinement applied to original mesh
     */
    jsoncons::optional<int> refineLevel{};
    /**
     * @brief Upper threshold for filtering (w.r.t mean)
     */
    jsoncons::optional<double> upperThreshold{};
    /**
     * @brief Lower threshold for filtering (w.r.t mean)
     */
    jsoncons::optional<double> lowerThreshold{};
    /**
     * @brief If @c periodic3D has value, then options for the 3D mesh.
     * Otherwise, only stl file is generated
     */
    jsoncons::optional<Periodic3DOpts> periodic3DOpts;

    JSONCONS_TYPE_TRAITS_FRIEND

   private:
    Opts() = default;
    static constexpr const char *type = "surface";
    static constexpr const char *engine = "packmesh";
  };

  SurfacePackMeshDriver(Files files, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files files);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  SurfacePackMeshDriver();

  Files files_;
  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_SURFACEPACKMESHDRIVER_H_
