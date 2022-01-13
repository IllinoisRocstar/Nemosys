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
#include "Drivers/MeshGen/SnappyMeshMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "Mesh/meshBase.H"
#include "MeshGeneration/snappymeshGen.H"

namespace NEM {
namespace DRV {

SnappyMeshMeshGenDriver::Opts::Opts(snappymeshParams params)
    : params(std::move(params)) {}

SnappyMeshMeshGenDriver::SnappyMeshMeshGenDriver(Files files,
                                                 snappymeshParams params)
    : files_(std::move(files)), opts_(std::move(params)) {}

SnappyMeshMeshGenDriver::SnappyMeshMeshGenDriver()
    : SnappyMeshMeshGenDriver({{}, {}}, {}) {}

const SnappyMeshMeshGenDriver::Files &SnappyMeshMeshGenDriver::getFiles()
    const {
  return files_;
}

void SnappyMeshMeshGenDriver::setFiles(Files files) {
  this->opts_.params.geomFileName = files.inputGeoFile;
  this->files_ = std::move(files);
}

const snappymeshParams &SnappyMeshMeshGenDriver::getParams() const {
  return getOpts().params;
}

void SnappyMeshMeshGenDriver::setParams(snappymeshParams params) {
  setOpts(Opts{std::move(params)});
}

const SnappyMeshMeshGenDriver::Opts &SnappyMeshMeshGenDriver::getOpts() const {
  return opts_;
}

void SnappyMeshMeshGenDriver::setOpts(Opts opts) {
  if (!opts.params.geomFileName.empty()) {
    this->files_.inputGeoFile = opts.params.geomFileName;
    this->opts_ = std::move(opts);
  } else {
    this->opts_ = std::move(opts);
    this->opts_.params.geomFileName = this->files_.inputGeoFile;
  }
}

void SnappyMeshMeshGenDriver::execute() const {
  auto paramsCopy = this->opts_.params;
  snappymeshGen generator{&paramsCopy};
  // TODO: Make sure blockMeshGen::createMeshFromSTL sets return value and check
  //  it here
  // Parameter not used
  generator.createMeshFromSTL(nullptr);
  std::string newname = nemAux::trim_fname(this->files_.inputGeoFile, ".vtu");
  auto mesh =
      meshBase::CreateShared(meshBase::Create(generator.getDataSet(), newname));
  mesh->setFileName(this->files_.outputMeshFile);
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
