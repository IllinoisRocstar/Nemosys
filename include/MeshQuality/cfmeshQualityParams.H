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
#ifndef NEMOSYS_CFMESHQUALITYPARAMS_H_
#define NEMOSYS_CFMESHQUALITYPARAMS_H_

#ifdef HAVE_CFMSH

#include "nemosys_export.h"
#include "MeshQuality/meshQualityParams.H"

#include <jsoncons/config/jsoncons_config.hpp>  // for jsoncons::optional
#include <jsoncons/json_traits_macros.hpp>
#include <string>

class NEMOSYS_EXPORT cfmeshQualityParams : public meshQualityParams {
 public:
  // initializing params with defaults
  cfmeshQualityParams() = default;

  int nIterations{50};          // number of optimization iterations
  int nLoops{10};               // number of inner loops in optimization
  double qualThrsh{0.1};         // quality threshold
  int nSrfItr{2};              // number of surface iterations
  jsoncons::optional<std::string> consCellSet{};  // name of constrained cellSet

  JSONCONS_TYPE_TRAITS_FRIEND

 protected:
  static constexpr const char *method = "meshOptimizer";
};

#endif  // HAVE_CFMSH

#endif  // NEMOSYS_CFMESHQUALITYPARAMS_H_
