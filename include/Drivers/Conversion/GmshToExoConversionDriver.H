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
#ifndef NEMOSYS_GMSHTOEXOCONVERSIONDRIVER_H_
#define NEMOSYS_GMSHTOEXOCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT GmshToExoConversionDriver : public ConversionDriver {
 public:
  using Files = DriverOutFile;

  struct NEMOSYS_EXPORT PostProcTask {
   public:
    explicit PostProcTask(std::string taskFile);
    std::string taskFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    PostProcTask() = default;
  };

  struct NEMOSYS_EXPORT MeshData {
   public:
    explicit MeshData(std::string meshFile);
    std::string meshFileName{};
    std::string meshName{"default"};
    bool usePhys{false};
    bool makeFreeSurfSS{false};
    bool splitTopBotSS{false};
    std::vector<std::string> sideSetNames{};
    std::map<std::string, std::string> elmBlkNames{};
    jsoncons::optional<std::string> addGlobalNodeSet{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    MeshData() = default;
  };

  struct NEMOSYS_EXPORT Opts {
   public:
    Opts(int numMeshes, std::vector<MeshData> meshData, bool needsPostProc);
    int numMeshes{};
    std::vector<MeshData> meshData{};
    bool needsPostProc{};
    int numTasks{0};
    std::vector<PostProcTask> tasks{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *methodName = "GMSH->EXO";
  };

  GmshToExoConversionDriver(Files file, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files file);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  GmshToExoConversionDriver();

  Files file_;
  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_GMSHTOEXOCONVERSIONDRIVER_H_
