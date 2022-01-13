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
#ifndef NEMOSYS_BLOCKMESHMESHGENDRIVER_H_
#define NEMOSYS_BLOCKMESHMESHGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/MeshGen/MeshGenDriver.H"
#include "MeshGeneration/blockMeshParams.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT BlockMeshMeshGenDriver : public MeshGenDriver {
 public:
  using Files = DriverOutFile;

  BlockMeshMeshGenDriver(Files file, blockMeshParams params);

  const Files &getFiles() const;
  void setFiles(Files file);
  const blockMeshParams &getParams() const;
  void setParams(blockMeshParams params);
  void execute() const override;

  struct NEMOSYS_NO_EXPORT Opts {
   public:
    explicit Opts(blockMeshParams params);
    blockMeshParams params{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *meshEngine = "blockMesh";
  };

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  BlockMeshMeshGenDriver();

  Files file_;
  Opts opts_;

  const Opts &getOpts() const;
  void setOpts(Opts opts);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_BLOCKMESHMESHGENDRIVER_H_
