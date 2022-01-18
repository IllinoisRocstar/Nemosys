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
#include "Drivers/MeshQualityDriver.H"

#ifdef HAVE_CFMSH
#  include "MeshQuality/MeshQuality.H"
#endif

namespace NEM {
namespace DRV {

jsoncons::string_view MeshQualityDriver::getProgramType() const {
  return programType;
}

CheckMeshQualDriver::Files::Files(std::string input, std::string output)
    : inputMeshFile(std::move(input)), outputFile(std::move(output)) {}

CheckMeshQualDriver::CheckMeshQualDriver(Files files)
    : files_(std::move(files)) {}

CheckMeshQualDriver::CheckMeshQualDriver() : CheckMeshQualDriver({{}, {}}) {}

CheckMeshQualDriver::Opts CheckMeshQualDriver::getOpts() { return {}; }

const CheckMeshQualDriver::Files &CheckMeshQualDriver::getFiles() const {
  return files_;
}

void CheckMeshQualDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

void CheckMeshQualDriver::execute() const {
  auto mesh = meshBase::Create(this->files_.inputMeshFile);
  mesh->checkMesh(this->files_.outputFile);
}

#ifdef HAVE_CFMSH
OptimizeMeshQualDriver::Opts::Opts(std::vector<cfmeshQualityParams> params)
    : params(std::move(params)) {}

OptimizeMeshQualDriver::OptimizeMeshQualDriver(
    std::vector<cfmeshQualityParams> params)
    : opts_(std::move(params)) {}

OptimizeMeshQualDriver::OptimizeMeshQualDriver()
    : OptimizeMeshQualDriver(std::vector<cfmeshQualityParams>{}) {}

const std::vector<cfmeshQualityParams> &OptimizeMeshQualDriver::getParams()
    const {
  return getOpts().params;
}

void OptimizeMeshQualDriver::setParams(std::vector<cfmeshQualityParams> params) {
  setOpts(Opts{std::move(params)});
}

void OptimizeMeshQualDriver::addParams(cfmeshQualityParams params) {
  this->opts_.params.emplace_back(std::move(params));
}

const OptimizeMeshQualDriver::Opts &OptimizeMeshQualDriver::getOpts() const {
  return opts_;
}

void OptimizeMeshQualDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void OptimizeMeshQualDriver::execute() const {
  for (const auto &param : this->opts_.params) {
    MeshQuality mq{&param};
    mq.cfmOptimize();
  }
}
#endif

}  // namespace DRV
}  // namespace NEM
