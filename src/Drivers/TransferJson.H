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
#ifndef NEMOSYS_TRANSFERJSON_H_
#define NEMOSYS_TRANSFERJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/TransferDriver.H"

// TransferDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::TransferDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, "Transfer Options", JSONCONS_RDWR),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::TransferDriver::programType;
     }))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::TransferDriver::Files,
                                (sourceMeshFile, "Source Mesh File"),
                                (targetMeshFile, "Target Mesh File"),
                                (outputMeshFile, NEM::DRV::JSON::outMeshFile))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::TransferDriver::Opts, 2,
                              (method, "Method"),
                              (checkQuality, "Check Transfer Quality"),
                              (arrayNames, "Array Names"))

#endif  // NEMOSYS_TRANSFERJSON_H_
