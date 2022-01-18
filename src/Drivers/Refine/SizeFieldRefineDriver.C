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
#include "Drivers/Refine/SizeFieldRefineDriver.H"

#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

SizeFieldRefineDriver::Opts::Opts(Method method, std::string arrayName,
                                  double stdDevMult, bool maxIsMin,
                                  bool transferData)
    : method(method),
      arrayName(std::move(arrayName)),
      stdDevMult(stdDevMult),
      maxIsMin(maxIsMin),
      transferData(transferData) {}

std::string SizeFieldRefineDriver::Opts::getMethodStr() const {
  switch (this->method) {
    case Method::VALUE: return valStr;
    case Method::GRADIENT: return gradStr;
  }
  return "";
}

SizeFieldRefineDriver::SizeFieldRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

SizeFieldRefineDriver::SizeFieldRefineDriver()
    : SizeFieldRefineDriver({{}, {}}, {{}, {}, {}, {}, {}}) {}

const SizeFieldRefineDriver::Opts &SizeFieldRefineDriver::getOpts() const {
  return opts_;
}

void SizeFieldRefineDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void SizeFieldRefineDriver::execute() const {
  std::cout << "Size Factor = " << this->opts_.sizeFactor << std::endl;
  std::shared_ptr<meshBase> mesh =
      meshBase::CreateShared(this->files_.inputMeshFile);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  // Edge scale is unused
  mesh->refineMesh(this->opts_.getMethodStr(), this->opts_.arrayName,
                   this->opts_.stdDevMult, this->opts_.maxIsMin, {},
                   this->files_.outputMeshFile, this->opts_.transferData,
                   this->opts_.sizeFactor);
}

}  // namespace DRV
}  // namespace NEM
