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
#include "Drivers/MeshGen/BlockMeshMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "MeshGeneration/blockMeshGen.H"
#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

BlockMeshMeshGenDriver::Opts::Opts(blockMeshParams params)
    : params(std::move(params)) {}

BlockMeshMeshGenDriver::BlockMeshMeshGenDriver(Files file,
                                               blockMeshParams params)
    : file_(std::move(file)), opts_(std::move(params)) {}

BlockMeshMeshGenDriver::BlockMeshMeshGenDriver()
    : BlockMeshMeshGenDriver(Files{{}}, {}) {}

const BlockMeshMeshGenDriver::Files &BlockMeshMeshGenDriver::getFiles() const {
  return file_;
}

void BlockMeshMeshGenDriver::setFiles(Files file) {
  this->file_ = std::move(file);
}

const blockMeshParams &BlockMeshMeshGenDriver::getParams() const {
  return getOpts().params;
}

void BlockMeshMeshGenDriver::setParams(blockMeshParams params) {
  setOpts(Opts{std::move(params)});
}

const BlockMeshMeshGenDriver::Opts &BlockMeshMeshGenDriver::getOpts() const {
  return opts_;
}

void BlockMeshMeshGenDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void BlockMeshMeshGenDriver::execute() const {
  auto paramsCopy = this->opts_.params;
  blockMeshGen generator{&paramsCopy};
  // TODO: Make sure blockMeshGen::createMeshFromSTL sets return value and check
  //  it here
  // Parameter not used
  generator.createMeshFromSTL(nullptr);
  auto mesh = meshBase::CreateShared(
      meshBase::Create(generator.getDataSet(), this->file_.outputFile));
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
