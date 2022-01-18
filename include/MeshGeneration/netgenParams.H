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
#ifndef NEMOSYS_NETGENPARAMS_H_
#define NEMOSYS_NETGENPARAMS_H_

#include "nemosys_export.h"
#include "MeshGeneration/meshingParams.H"

#include <string>

class NEMOSYS_EXPORT netgenParams : public meshingParams
{
  public:
    // initializing params with defaults
    netgenParams() = default;

  public:
    bool uselocalh{true};
    double maxh{1000.0};
    double fineness{0.5};
    double grading{0.3};
    double elementsperedge{2.0};
    double elementspercurve{2.0};
    bool closeedgeenable{false};
    double closeedgefact{2.0};
    bool second_order{false};
    std::string meshsize_filename{};
    bool quad_dominated{false};
    bool optvolmeshenable{true};
    int optsteps_2d{3};
    int optsteps_3d{3};
    bool invert_tets{false};
    bool invert_trigs{false};
    bool check_overlap{true};
    bool check_overlapping_boundary{true};
    bool refine_with_geom{false};
    bool refine_without_geom{false};
};

#endif  // NEMOSYS_NETGENPARAMS_H_
