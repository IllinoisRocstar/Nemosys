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
#ifndef NEMOSYS_MESHQUALITYJSON_H_
#define NEMOSYS_MESHQUALITYJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/MeshQualityDriver.H"

#ifdef HAVE_CFMSH
#  define IF_CFMSH(a) a,
#else
#  define IF_CFMSH(a)
#endif

namespace NEM {
namespace DRV {
namespace JSON {

static constexpr auto meshQualOpts = "Mesh Quality Options";
static constexpr auto meshQualEngine = "Mesh Quality Engine";

}  // namespace JSON
}  // namespace DRV
}  // namespace NEM

// MeshQualityDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_INTERMEDIATE(
    NEM::DRV::MeshQualityDriver,
    (IF_CFMSH(NEM::DRV::OptimizeMeshQualDriver)
     NEM::DRV::CheckMeshQualDriver),
    NEM::DRV::NemDriver, 1,
    (, , NEM::DRV::JSON::meshQualOpts, NEM_JSON_RDONLY_OVERRIDE,
     NEM_JSON_CHECK_KEY_ONLY),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::MeshQualityDriver::programType;
     }))

// CheckMeshQualDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::CheckMeshQualDriver, NEM::DRV::MeshQualityDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, , NEM::DRV::JSON::meshQualOpts, JSONCONS_RDONLY))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::CheckMeshQualDriver::Files,
                                (inputMeshFile, NEM::DRV::JSON::inMeshFile),
                                (outputFile, "Output File"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::CheckMeshQualDriver::Opts,
    (meshQualityEngine, NEM::DRV::JSON::meshQualEngine, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::CheckMeshQualDriver::Opts::meshQualityEngine;
     }))

#ifdef HAVE_CFMSH
// OptimizeMeshQualDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::OptimizeMeshQualDriver, NEM::DRV::MeshQualityDriver, 1,
    (getOpts, setOpts, NEM::DRV::JSON::meshQualOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::OptimizeMeshQualDriver::Opts,
    (meshQualityEngine, NEM::DRV::JSON::meshQualEngine, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::OptimizeMeshQualDriver::Opts::meshQualityEngine;
     }),
    (params, "Schedule"))

JSONCONS_N_MEMBER_NAME_TRAITS(cfmeshQualityParams, 5,
                              (method, "Method", JSONCONS_RDONLY,
                               [](const jsoncons::string_view &x) {
                                 return x == cfmeshQualityParams::method;
                               }),
                              (nIterations, "NIterations"), (nLoops, "NLoops"),
                              (qualThrsh, "QualityThreshold"),
                              (nSrfItr, "NSurfaceIterations"),
                              (consCellSet, "ConstrainedCellSet"))
#endif

#undef IF_CFMSH

#endif  // NEMOSYS_MESHQUALITYJSON_H_
