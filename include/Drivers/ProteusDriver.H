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
#ifndef NEMOSYS_PROTEUSDRIVER_H_
#define NEMOSYS_PROTEUSDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

namespace NEM {
namespace DRV {

/**
 * @brief Driver to read Proteus output. Converts Proteus HDF mesh to VTK/Exodus
 * format
 */
class NEMOSYS_EXPORT ProteusDriver : public NemDriver {
 public:
  struct NEMOSYS_EXPORT Files {
   public:
    Files(std::string fieldFName, std::string meshFName,
          std::string exoMeshFName);
    /**
     * @brief Proteus format HDF5 field filename
     */
    std::string fieldFName{};
    /**
     * @brief Output VTK mesh filename
     */
    std::string meshFName{};
    /**
     * @brief Output Exodus mesh filename
     */
    std::string exoMeshFName{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Files() = default;
  };

  struct NEMOSYS_EXPORT Opts {
   public:
    explicit Opts(std::string edgeSidesetName);
    /**
     * @brief Name of sideset written to output Exodus file
     */
    std::string edgeSidesetName{};
    /**
     * @brief Boolean converting high order cells to low order
     */
    bool lowOrder{false};
    /**
     * @brief Boolean to employ boundary constraint during refinement
     */
    bool bndryConst{true};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };
  ProteusDriver(Files files, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files files);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  ~ProteusDriver() override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  ProteusDriver();

  Files files_;
  Opts opts_;

  static constexpr const char *programType = "Proteus";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_PROTEUSDRIVER_H_
