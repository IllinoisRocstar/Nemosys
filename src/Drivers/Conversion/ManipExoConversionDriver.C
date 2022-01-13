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
#include "Drivers/Conversion/ManipExoConversionDriver.H"

#include "Mesh/exoMesh.H"

namespace NEM {
namespace DRV {

ManipExoConversionDriver::ManipExoConversionDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

ManipExoConversionDriver::ManipExoConversionDriver()
    : ManipExoConversionDriver(Files{{}, {}}, Opts{{}}) {}

const ManipExoConversionDriver::Files &ManipExoConversionDriver::getFiles()
    const {
  return files_;
}

void ManipExoConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const ManipExoConversionDriver::Opts &ManipExoConversionDriver::getOpts()
    const {
  return opts_;
}

void ManipExoConversionDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void ManipExoConversionDriver::execute() const {
  // read in exo mesh from file name
  NEM::MSH::EXOMesh::exoMesh em{this->files_.inputMeshFile};
  em.read(this->files_.inputMeshFile);

  // Combining Blocks
  for (const auto &cmb : this->opts_.combineBlocks) {
    em.combineElmBlks(cmb.blockIds, cmb.newName);
  }

  // Set the output file name
  em.setFileName(this->files_.outputMeshFile);
  em.write();
  em.report();
}

}  // namespace DRV
}  // namespace NEM
