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
#include "Drivers/ProteusDriver.H"

#include "Mesh/meshBase.H"
#include "IO/proteusHdf5.H"

namespace NEM {
namespace DRV {

ProteusDriver::Files::Files(std::string fieldFName, std::string meshFName,
                            std::string exoMeshFName)
    : fieldFName(std::move(fieldFName)),
      meshFName(std::move(meshFName)),
      exoMeshFName(std::move(exoMeshFName)) {}

ProteusDriver::Opts::Opts(std::string edgeSidesetName)
    : edgeSidesetName(std::move(edgeSidesetName)) {}

ProteusDriver::ProteusDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {
  std::cout << "ProteusDriver created\n";
}

ProteusDriver::ProteusDriver() : ProteusDriver({{}, {}, {}}, Opts{{}}) {}

const ProteusDriver::Files &ProteusDriver::getFiles() const { return files_; }

void ProteusDriver::setFiles(Files files) { this->files_ = std::move(files); }

const ProteusDriver::Opts &ProteusDriver::getOpts() const { return opts_; }

void ProteusDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

ProteusDriver::~ProteusDriver() {
  std::cout << "ProteusDriver destroyed" << std::endl;
}

jsoncons::string_view ProteusDriver::getProgramType() const {
  return programType;
}

void ProteusDriver::execute() const {
  proteusHdf5(this->files_.fieldFName, this->files_.meshFName,
              this->opts_.edgeSidesetName, this->files_.exoMeshFName,
              this->opts_.lowOrder, this->opts_.bndryConst);
}

}  // namespace DRV
}  // namespace NEM
