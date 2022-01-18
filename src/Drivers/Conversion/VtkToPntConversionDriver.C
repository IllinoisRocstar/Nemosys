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
#include "Drivers/Conversion/VtkToPntConversionDriver.H"

namespace NEM {
namespace DRV {

VtkToPntConversionDriver::Opts::Opts(int dim, PNTMesh::BlockMap blockMap)
    : dim(dim), elemBlockMap(std::move(blockMap)) {}

VtkToPntConversionDriver::VtkToPntConversionDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

VtkToPntConversionDriver::VtkToPntConversionDriver()
    : VtkToPntConversionDriver({{}, {}}, {{}, {}}) {}

const VtkToPntConversionDriver::Files &VtkToPntConversionDriver::getFiles()
    const {
  return files_;
}

void VtkToPntConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const VtkToPntConversionDriver::Opts &VtkToPntConversionDriver::getOpts()
    const {
  return opts_;
}

void VtkToPntConversionDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void VtkToPntConversionDriver::execute() const {
  auto source = meshBase::Create(files_.inputMeshFile);
  std::cout << "Number of Blocks : " << opts_.elemBlockMap.size() << std::endl;
  auto *pm = new PNTMesh::pntMesh(source, opts_.dim, opts_.elemBlockMap.size(),
                                  opts_.elemBlockMap);
  pm->write(files_.outputMeshFile);
  delete pm;
}

}  // namespace DRV
}  // namespace NEM
