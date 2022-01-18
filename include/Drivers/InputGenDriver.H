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
#ifndef NEMOSYS_INPUTGENDRIVER_H_
#define NEMOSYS_INPUTGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

#include <string>
#include <jsoncons/json.hpp>

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT InputGenDriver : public NemDriver {
 public:
  InputGenDriver(std::string service, jsoncons::json opts);

  const std::string &getService() const;
  void setService(std::string service);
  const jsoncons::json &getOpts() const;
  void setOpts(jsoncons::json opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  InputGenDriver() = default;

  static constexpr const char *programType = "Input Generation";

  std::string service_{};
  jsoncons::json opts_{};

  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_INPUTGENDRIVER_H_
