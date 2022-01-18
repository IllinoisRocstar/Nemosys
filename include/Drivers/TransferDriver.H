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
#ifndef NEMOSYS_TRANSFERDRIVER_H_
#define NEMOSYS_TRANSFERDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"
#include "Transfer/TransferBase.H"
#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT TransferDriver : public NemDriver {
 public:
  struct NEMOSYS_EXPORT Files {
   public:
    Files(std::string source, std::string target, std::string output);
    std::string sourceMeshFile{};
    std::string targetMeshFile{};
    std::string outputMeshFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Files() = default;
  };

  struct NEMOSYS_EXPORT Opts {
   public:
    Opts(std::string method, bool checkQuality);
    std::string method{};
    bool checkQuality{};
    jsoncons::optional<std::vector<std::string>> arrayNames{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };

  TransferDriver(Files files, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files files);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  ~TransferDriver() override;

  static std::shared_ptr<TransferBase> CreateTransferObject(
      meshBase *srcmsh, meshBase *trgmsh, const std::string &method);

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  TransferDriver();

  static constexpr const char *programType = "Transfer";
  jsoncons::string_view getProgramType() const override;

  Files files_;
  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_TRANSFERDRIVER_H_
