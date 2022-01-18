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
#ifndef NEMOSYS_PROTEUSJSON_H_
#define NEMOSYS_PROTEUSJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/ProteusDriver.H"

// ProteusDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::ProteusDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, "Proteus Options", JSONCONS_RDWR),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::ProteusDriver::programType;
     }))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::ProteusDriver::Files,
                                (fieldFName, "HDF5 Field File"),
                                (meshFName, NEM::DRV::JSON::outMeshFile),
                                (exoMeshFName, "Output Exodus File"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::ProteusDriver::Opts, 1,
                              (edgeSidesetName, "Edge Sideset"),
                              (lowOrder, "Output Low Order"),
                              (bndryConst, "Constraint Boundary"))

#endif  // NEMOSYS_PROTEUSJSON_H_
