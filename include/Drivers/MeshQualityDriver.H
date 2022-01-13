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
#ifndef NEMOSYS_MESHQUALITYDRIVER_H_
#define NEMOSYS_MESHQUALITYDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

#include "Mesh/meshBase.H"
#ifdef HAVE_CFMSH
#  include "MeshQuality/cfmeshQualityParams.H"
#endif

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT MeshQualityDriver : public NemDriver {
  JSONCONS_TYPE_TRAITS_FRIEND
 private:
  static constexpr const char *programType = "Mesh Quality";
  jsoncons::string_view getProgramType() const override;
};

class NEMOSYS_EXPORT CheckMeshQualDriver : public MeshQualityDriver {
 public:
  struct NEMOSYS_EXPORT Files {
   public:
    Files(std::string input, std::string output);
    std::string inputMeshFile{};
    std::string outputFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Files() = default;
  };
  explicit CheckMeshQualDriver(Files files);

  const Files &getFiles() const;
  void setFiles(Files files);

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  CheckMeshQualDriver();

  Files files_;

  void execute() const override;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *meshQualityEngine = "default";
  };

 private:
  static Opts getOpts();
};

#ifdef HAVE_CFMSH

class NEMOSYS_EXPORT OptimizeMeshQualDriver : public MeshQualityDriver {
 public:
  explicit OptimizeMeshQualDriver(std::vector<cfmeshQualityParams> params);

  const std::vector<cfmeshQualityParams> &getParams() const;
  void setParams(std::vector<cfmeshQualityParams> params);
  void addParams(cfmeshQualityParams params);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  OptimizeMeshQualDriver();

 public:
  struct NEMOSYS_NO_EXPORT Opts {
   public:
    explicit Opts(std::vector<cfmeshQualityParams> params);
    // Change this to std::vector<std::shared_ptr<meshQualityParams>> (or
    // std::vector<std::variant<>>) to support multiple methods in the future
    std::vector<cfmeshQualityParams> params{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *meshQualityEngine = "cfmesh";
  };

 private:
  Opts opts_;
  const Opts &getOpts() const;
  void setOpts(Opts opts);
};
#endif

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_MESHQUALITYDRIVER_H_
