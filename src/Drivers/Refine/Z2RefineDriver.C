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
#include "Drivers/Refine/Z2RefineDriver.H"

#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

Z2RefineDriver::Opts::Opts(bool transferData, std::string arrayName, int order)
    : transferData(transferData),
      arrayName(std::move(arrayName)),
      order(order) {}

Z2RefineDriver::Z2RefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

Z2RefineDriver::Z2RefineDriver() : Z2RefineDriver({{}, {}}, {{}, {}, {}}) {}

const Z2RefineDriver::Opts &Z2RefineDriver::getOpts() const {
  return opts_;
}

void Z2RefineDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void Z2RefineDriver::execute() const {
  std::shared_ptr<meshBase> mesh =
      meshBase::CreateShared(this->files_.inputMeshFile);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(Opts::method, this->opts_.arrayName, this->opts_.order,
                   this->files_.outputMeshFile, this->opts_.transferData);
}

}  // namespace DRV
}  // namespace NEM
