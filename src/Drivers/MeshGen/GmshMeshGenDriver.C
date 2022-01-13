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
#include "Drivers/MeshGen/GmshMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "MeshGeneration/gmshGen.H"
#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

GmshMeshGenDriver::Opts::Opts(NEM::GEN::gmshParams params)
    : params(std::move(params)) {}

GmshMeshGenDriver::GmshMeshGenDriver(Files files, NEM::GEN::gmshParams params)
    : files_(std::move(files)), opts_(std::move(params)) {}

GmshMeshGenDriver::GmshMeshGenDriver() : GmshMeshGenDriver({{}, {}}, {}) {}

const GmshMeshGenDriver::Files &GmshMeshGenDriver::getFiles() const {
  return files_;
}

void GmshMeshGenDriver::setFiles(Files files) {
  this->opts_.params.ofname = files.outputMeshFile;
  this->files_ = std::move(files);
}

const NEM::GEN::gmshParams &GmshMeshGenDriver::getParams() const {
  return getOpts().params;
}

void GmshMeshGenDriver::setParams(NEM::GEN::gmshParams params) {
  setOpts(Opts{std::move(params)});
}

const GmshMeshGenDriver::Opts &GmshMeshGenDriver::getOpts() const {
  return opts_;
}

void GmshMeshGenDriver::setOpts(Opts opts) {
  if (!opts.params.ofname.empty()) {
    this->files_.outputMeshFile = opts.params.ofname;
    this->opts_ = std::move(opts);
  } else {
    this->opts_ = std::move(opts);
    this->opts_.params.ofname = this->files_.outputMeshFile;
  }
}

void GmshMeshGenDriver::execute() const {
  auto paramsCopy = this->opts_.params;
  NEM::GEN::gmshGen generator{&paramsCopy};
  int status = generator.createMeshFromSTL(this->files_.inputGeoFile.c_str());
  if (status) {
    std::cerr << "Mesh Generation encountered error." << std::endl;
    exit(1);
  }
  std::string outputType = nemAux::find_ext(this->files_.outputMeshFile);
  if (outputType == ".msh") {
    return;
  }
  std::string newname = nemAux::trim_fname(this->files_.inputGeoFile, ".msh");
  auto mesh = meshBase::CreateShared(meshBase::exportGmshToVtk(newname));
  mesh->setFileName(this->files_.outputMeshFile);
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
