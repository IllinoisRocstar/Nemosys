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
#ifndef NEMOSYS_Z2REFINEDRIVER_H_
#define NEMOSYS_Z2REFINEDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Refine/RefineDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT Z2RefineDriver : public RefineDriver {
 public:
  struct NEMOSYS_EXPORT Opts {
   public:
    Opts(bool transferData, std::string arrayName, int order);
    bool transferData{};
    std::string arrayName{};
    int order{};
    static constexpr const char *method = "Z2 Error Estimator";
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };
  Z2RefineDriver(Files files, Opts opts);

  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  Z2RefineDriver();

  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_Z2REFINEDRIVER_H_
