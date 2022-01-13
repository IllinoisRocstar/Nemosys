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
#include "Drivers/MeshGen/NetgenMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "Mesh/meshBase.H"
#include "MeshGeneration/netgenGen.H"

namespace NEM {
namespace DRV {

NetgenMeshGenDriver::Opts::Opts(netgenParams params)
    : params(std::move(params)) {}

NetgenMeshGenDriver::NetgenMeshGenDriver(Files files, netgenParams params)
    : files_(std::move(files)), opts_(std::move(params)) {}

NetgenMeshGenDriver::NetgenMeshGenDriver()
    : NetgenMeshGenDriver({{}, {}}, {}) {}

const NetgenMeshGenDriver::Files &NetgenMeshGenDriver::getFiles() const {
  return files_;
}

void NetgenMeshGenDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const netgenParams &NetgenMeshGenDriver::getParams() const {
  return getOpts().params;
}

void NetgenMeshGenDriver::setParams(netgenParams params) {
  setOpts(Opts{std::move(params)});
}

const NetgenMeshGenDriver::Opts &NetgenMeshGenDriver::getOpts() const {
  return opts_;
}

void NetgenMeshGenDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

void NetgenMeshGenDriver::execute() const {
  netgenGen generator{&this->opts_.params};
  int status = generator.createMeshFromSTL(this->files_.inputGeoFile.c_str());
  if (status) {
    std::cerr << "Mesh Generation encountered error." << std::endl;
    exit(1);
  }
  std::string newname = nemAux::trim_fname(this->files_.inputGeoFile, ".vol");
  auto mesh = meshBase::CreateShared(meshBase::exportVolToVtk(newname));
  mesh->setFileName(this->files_.outputMeshFile);
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
