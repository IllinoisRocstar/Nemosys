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
#ifndef NEMOSYS_AUTOVERIFICATIONDRIVER_H_
#define NEMOSYS_AUTOVERIFICATIONDRIVER_H_

#include "nemosys_export.h"

#include "Drivers/NemDriver.H"

namespace NEM {
namespace DRV {

/**
 * This class drives the automatic verification (AV) in Nemosys. Three
 * successively refined meshes of a given geometry are required to evaluate
 * whether the coarsest mesh is within the asymptotic range. The grid
 * convergence index (GCI) of each component of each selected field is
 * evaluated: if they are all approximately equal to 1 then the coarse mesh is
 * sufficiently refined for an accurate solution (and is within asymptotic
 * range).
 */
class NEMOSYS_EXPORT AutoVerificationDriver : public NemDriver {
 public:
  struct NEMOSYS_EXPORT Files {
   public:
    /**
     * @brief Set AV driver files
     * @param coarseMesh coarsest mesh
     * @param fineMesh fine mesh
     * @param finerMesh finest mesh
     */
    Files(std::string coarseMesh, std::string fineMesh,
          std::string finerMesh);
    std::string coarseMeshFile{};
    std::string fineMeshFile{};
    std::string finerMeshFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Files() = default;
  };

  struct NEMOSYS_EXPORT Opts {
   public:
    /**
     * @brief Set AV driver options
     * @param arrayIDs Set @c arrayIds member
     */
    explicit Opts(std::vector<int> arrayIds);
    /**
     * @brief field ids for which we will evaluate the grid convergence index
     * (GCI)
     */
    std::vector<int> arrayIds;
    /**
     * @brief default is 'Consistent Interpolation',
     * alternatives include 'Conservative Surface Transfer' (if enabled,
     * requires IMPACT) and 'Conservative Volume Transfer' (if enabled)
     */
    std::string transferType{"Consistent Interpolation"};
    /**
     * @brief the target GCI value, set to 1.1 by default
     */
    double targetGCI{1.1};
    /**
     * @brief sets number of threads used in transfer (if OpenMP enabled); if
     * unset and OpenMP enabled, @c omp_get_max_threads used
     */
    jsoncons::optional<int> numThreads{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };

  AutoVerificationDriver(Files files, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files files);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  AutoVerificationDriver();
  static constexpr const char *programType = "Verification";
  jsoncons::string_view getProgramType() const override;

  Files files_;
  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_AUTOVERIFICATIONDRIVER_H_
