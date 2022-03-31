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
#ifndef NEMOSYS_CONVERSIONJSON_H_
#define NEMOSYS_CONVERSIONJSON_H_

#include <jsoncons/json.hpp>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/Conversion/ConversionDriver.H"
#include "Drivers/Conversion/GmshToExoConversionDriver.H"
#include "Drivers/Conversion/GmshToVtkConversionDriver.H"
#include "Drivers/Conversion/ManipExoConversionDriver.H"
#include "Drivers/Conversion/SmartConversionDriver.H"
#include "Drivers/Conversion/VtkHexToTetConversionDriver.H"
#include "Drivers/Conversion/VtkToCobaltConversionDriver.H"
#include "Drivers/Conversion/VtkToPatranConversionDriver.H"
#include "Drivers/Conversion/VtkToPntConversionDriver.H"

#ifdef HAVE_CFMSH
#  define IF_CFMSH(a) a,
#  include "Drivers/Conversion/FoamToVtkConversionDriver.H"
#  include "Drivers/Conversion/VtkToFoamConversionDriver.H"
#else
#  define IF_CFMSH(a)
#endif

#if defined(HAVE_CFMSH) && defined(HAVE_GMSH)
#  define IF_CFMSH_GMSH(a) a,
#  include "Drivers/Conversion/FoamToMshConversionDriver.H"
#else
#  define IF_CFMSH_GMSH(a)
#endif

namespace NEM {
namespace DRV {
namespace JSON {

static constexpr auto convOptions = "Conversion Options";
static constexpr auto convMethod = "Method";

inline std::vector<PNTMesh::surfaceBCTag> pntMeshBcTagFromJSON(
    const std::string &input) {
  return {PNTMesh::bcTagNum(input)};
}

inline std::string pntMeshBcTagToJSON(
    const std::vector<PNTMesh::surfaceBCTag> &input) {
  return input.empty() ? std::string{} : PNTMesh::bcTagStr(input.at(0));
}

inline std::vector<int> rangeToVec(const std::pair<int, int> &input) {
  std::vector<int> out(input.second - input.first + 1);
  std::iota(out.begin(), out.end(), input.first);
  return out;
}

}  // namespace JSON
}  // namespace DRV
}  // namespace NEM

// ConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_INTERMEDIATE(
    NEM::DRV::ConversionDriver,
    (IF_CFMSH_GMSH(NEM::DRV::FoamToMshConversionDriver)
     IF_CFMSH(NEM::DRV::FoamToVtkConversionDriver)
     IF_CFMSH(NEM::DRV::VtkToFoamConversionDriver)
     NEM::DRV::GmshToExoConversionDriver, NEM::DRV::GmshToVtkConversionDriver,
     NEM::DRV::ManipExoConversionDriver, NEM::DRV::SmartConversionDriver,
     NEM::DRV::VtkHexToTetConversionDriver,
     NEM::DRV::VtkToCobaltConversionDriver,
     NEM::DRV::VtkToPatranConversionDriver, NEM::DRV::VtkToPntConversionDriver),
    NEM::DRV::NemDriver, 1,
    (, , NEM::DRV::JSON::convOptions, NEM_JSON_RDONLY_OVERRIDE,
     NEM_JSON_CHECK_KEY_ONLY),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::ConversionDriver::programType;
     }))

#ifdef HAVE_CFMSH
#ifdef HAVE_GMSH
// FoamToMshConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::FoamToMshConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, , NEM::DRV::JSON::convOptions, JSONCONS_RDONLY))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::FoamToMshConversionDriver::Opts,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::FoamToMshConversionDriver::Opts::methodName;
     }))
#endif

// FoamToVtkConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::FoamToVtkConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, , NEM::DRV::JSON::convOptions, JSONCONS_RDONLY))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::FoamToVtkConversionDriver::Opts,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::FoamToVtkConversionDriver::Opts::methodName;
     }))

// VtkToFoamConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::VtkToFoamConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, , NEM::DRV::JSON::convOptions, JSONCONS_RDONLY))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::VtkToFoamConversionDriver::Files,
                                (inputMeshFile, NEM::DRV::JSON::inMeshFile))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::VtkToFoamConversionDriver::Opts,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::VtkToFoamConversionDriver::Opts::methodName;
     }))
#endif

// GmshToExoConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::GmshToExoConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, NEM::DRV::JSON::convOptions, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::GmshToExoConversionDriver::Opts, 4,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::GmshToExoConversionDriver::Opts::methodName;
     }),
    (numMeshes, "Number of Mesh"), (meshData, "Mesh Data"),
    (needsPostProc, "Post Processing"), (numTasks, "Number of Tasks"),
    (tasks, "Tasks"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::GmshToExoConversionDriver::MeshData, 1,
    (meshFileName, "File"), (meshName, "Name"),
    (usePhys, "Use Physical Groups"),
    (makeFreeSurfSS, "Free Surface Side Set"),
    (splitTopBotSS, "Split Top and Bottom"),
    (sideSetNames, "Side Set Names"),
    (elmBlkNames, "Element Block Names"),
    (addGlobalNodeSet, "Add Global Node Set"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::GmshToExoConversionDriver::PostProcTask, (taskFile, "File"))

// GmshToVtkConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::GmshToVtkConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, , NEM::DRV::JSON::convOptions, JSONCONS_RDONLY))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::GmshToVtkConversionDriver::Opts,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::GmshToVtkConversionDriver::Opts::methodName;
     }))

// ManipExoConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::ManipExoConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, NEM::DRV::JSON::convOptions, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::ManipExoConversionDriver::Opts, 1,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::ManipExoConversionDriver::Opts::methodName;
     }),
    (combineBlocks, "CombineBlocks"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::ManipExoConversionDriver::CombineBlocks,
                              0, (blockIds, "BlockIDs"), (newName, "Rename"))

// SmartConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::SmartConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, NEM::DRV::JSON::convOptions, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::SmartConversionDriver::Opts, 1,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::SmartConversionDriver::Opts::methodName;
     }),
    (meshData, "Mesh Data"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::SmartConversionDriver::MeshData, 0,
    (inputFormat, "Input Mesh Format"),
    (outputFormat, "Output Mesh Format"))

// VtkHexToTetConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::VtkHexToTetConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, , NEM::DRV::JSON::convOptions, JSONCONS_RDONLY))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::VtkHexToTetConversionDriver::Opts,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::VtkHexToTetConversionDriver::Opts::methodName;
     }))

// VtkToCobaltConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::VtkToCobaltConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, , NEM::DRV::JSON::convOptions, JSONCONS_RDONLY))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::VtkToCobaltConversionDriver::Files,
                                (inputMeshFile, NEM::DRV::JSON::inMeshFile),
                                (outputCgrFile, NEM::DRV::JSON::outMeshFile),
                                (outputCgiFile, "Output Patch Map File"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::VtkToCobaltConversionDriver::Opts,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::VtkToCobaltConversionDriver::Opts::methodName;
     }))

// VtkToPatranConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::VtkToPatranConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, NEM::DRV::JSON::convOptions, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::VtkToPatranConversionDriver::Opts, 3,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::VtkToPatranConversionDriver::Opts::methodName;
     }),
    (bcInfo, "BC Info"), (nodePatchPreference, "Node Patch Preference"))

NEM_JSON_N_MEMBER_NAME_TRAITS_BASE(
    NEM::DRV::VtkToPatranConversionDriver::BoundaryCond,
    (NEM::DRV::VtkToPatranConversionDriver::FaceBC,
     NEM::DRV::VtkToPatranConversionDriver::NodeBC),
    1, (patchNum, "Patch Number"), (getBCType(), "BC Type", JSONCONS_RDONLY))

NEM_JSON_WRAP_SMART_PTR(NEM::DRV::VtkToPatranConversionDriver::BoundaryCond,
                        std::shared_ptr)

NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(
    NEM::DRV::VtkToPatranConversionDriver::FaceBC,
    NEM::DRV::VtkToPatranConversionDriver::BoundaryCond, 1,
    (rocFracFSIType, "Rocfrac FSI Type"),
    (bcType, "BC Type", NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::VtkToPatranConversionDriver::FaceBC::bcType;
     }))

NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(
    NEM::DRV::VtkToPatranConversionDriver::NodeBC,
    NEM::DRV::VtkToPatranConversionDriver::BoundaryCond, 4,
    (rocfracControlType, "RocfracControl Type"), (structural, "Structural"),
    (meshMotion, "Mesh Motion"), (thermal, "Thermal"),
    (bcType, "BC Type", NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::VtkToPatranConversionDriver::NodeBC::bcType;
     }))

// VtkToPntConversionDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::VtkToPntConversionDriver, NEM::DRV::ConversionDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, NEM::DRV::JSON::convOptions, JSONCONS_RDWR))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::VtkToPntConversionDriver::Opts, 3,
    (methodName, NEM::DRV::JSON::convMethod, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::VtkToPntConversionDriver::Opts::methodName;
     }),
    (dim, "Dimension"), (elemBlockMap, "Block Data"))

// Function-like macro so that the Into parameter of the first elmIds clause
// knows the type to deserialize and then pass to rangeToVec. Not that the Into
// function is unused due to NEM_JSON_RDWR_OVERRIDE - so that the second elmIds
// "Element ID" clause is used to serialize
#define PAIR_INT_DUMMY_FUNC(X) \
  std::pair<int, int> {}
NEM_JSON_N_MEMBER_NAME_TRAITS_VAL(
    PNTMesh::blockType, 5, (ordIntrp, "Element Order"),
    (ordEquat, "Equation Order"),
    (eTpe, "Element Type", JSONCONS_RDWR, NEM_JSON_CHECK_KEY_ONLY,
     PNTMesh::elmTypeStr, PNTMesh::elmTypeNum),
    (srfBCTag, "BC Tag", JSONCONS_RDWR, NEM_JSON_CHECK_KEY_ONLY,
     NEM::DRV::JSON::pntMeshBcTagToJSON, NEM::DRV::JSON::pntMeshBcTagFromJSON),
    (regionName, "Name"),
    (elmIds, "Element ID Range", NEM_JSON_RDWR_OVERRIDE,
     NEM_JSON_CHECK_KEY_ONLY, PAIR_INT_DUMMY_FUNC, NEM::DRV::JSON::rangeToVec),
    (elmIds, "Element ID"))
#undef PAIR_INT_DUMMY_FUNC

#undef IF_CFMSH
#undef IF_CFMSH_GMSH

#endif  // NEMOSYS_CONVERSIONJSON_H_
