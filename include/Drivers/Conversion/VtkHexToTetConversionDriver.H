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
#ifndef NEMOSYS_VTKHEXTOTETCONVERSIONDRIVER_H_
#define NEMOSYS_VTKHEXTOTETCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT VtkHexToTetConversionDriver : public ConversionDriver {
 public:
  using Files = DriverInOutFiles;

  explicit VtkHexToTetConversionDriver(Files files);

  const Files &getFiles() const;
  void setFiles(Files files);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  VtkHexToTetConversionDriver();

  Files files_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "VTK-HEX->VTK-TET";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_VTKHEXTOTETCONVERSIONDRIVER_H_
