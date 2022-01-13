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
#include "Drivers/InputGenDriver.H"

#include <iostream>
#include "AuxiliaryFunctions.H"
#include "InputGeneration/ep16Post.H"
#include "InputGeneration/ep16Prep.H"

namespace NEM {
namespace DRV {

InputGenDriver::InputGenDriver(std::string service, jsoncons::json opts)
    : service_(std::move(service)), opts_(std::move(opts)) {}

const std::string &InputGenDriver::getService() const { return service_; }

void InputGenDriver::setService(std::string service) {
  this->service_ = std::move(service);
}

const jsoncons::json &InputGenDriver::getOpts() const { return opts_; }

void InputGenDriver::setOpts(jsoncons::json opts) {
  this->opts_ = std::move(opts);
}

jsoncons::string_view InputGenDriver::getProgramType() const {
  return programType;
}

void InputGenDriver::execute() const {
  std::string srvName = this->service_;
  nemAux::toLower(srvName);
  nemAux::Timer T;
  T.start();
  if (srvName == "epic_2016") {
#ifdef HAVE_EPIC
    ep16Prep::readJSON(this->opts_);
#else
    std::cerr << "Compile the code with EPIC module enabled." << std::endl;
    exit(0);
#endif
  } else if (srvName == "epic_2016_post") {
#ifdef HAVE_EPIC
    int ret;
    NEM::EPC::ep16Post::readJSON(this->opts_, ret);
#else
    std::cerr << "Compile the code with EPIC module enabled." << std::endl;
    exit(0);
#endif
  } else {
    std::cerr << "The input generation service " << srvName << "is unsupported."
              << std::endl;
  }
  T.stop();
}

}  // namespace DRV
}  // namespace NEM
