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
#ifndef NEMOSYS_DRIVERJSONTYPETRAITS_H_
#define NEMOSYS_DRIVERJSONTYPETRAITS_H_

#include <jsoncons/json.hpp>
#include <memory>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/AutoVerificationDriver.H"
#include "Drivers/Conversion/ConversionDriver.H"
#include "Drivers/MeshGen/MeshGenDriver.H"
#include "Drivers/MeshQualityDriver.H"
#include "Drivers/NemDriver.H"
#include "Drivers/PackMesh/PackMeshDriver.H"
#include "Drivers/Refine/RefineDriver.H"
#include "Drivers/TransferDriver.H"

#ifdef HAVE_CGNS
#  define IF_CGNS(a) a,
#  include "Drivers/RocPartCommGenDriver.H"
#else
#  define IF_CGNS(a)
#endif

#ifdef HAVE_EPIC
#  define IF_EPIC(a) a,
#  include "Drivers/InputGenDriver.H"
#else
#  define IF_EPIC(a)
#endif

#if defined(HAVE_HDF5) && defined(HAVE_GMSH)
#  define IF_HDF5_GMSH(a) a,
#  include "Drivers/ProteusDriver.H"
#else
#  define IF_HDF5_GMSH(a)
#endif

#ifdef HAVE_TEMPLATE_MESH
#  define IF_TEMPLATE_MESH(a) a,
#  include "Drivers/TemplateMeshDriver.H"
#else
#  define IF_TEMPLATE_MESH(a)
#endif

#ifdef HAVE_OCC
#  define IF_OCC(a) a,
#  include "Drivers/NucMeshDriver.H"
#else
#  define IF_OCC(a)
#endif

namespace NEM {
namespace DRV {
namespace JSON {

static constexpr auto programType = "Program Type";
static constexpr auto meshFiles = "Mesh File Options";
static constexpr auto inMeshFile = "Input Mesh File";
static constexpr auto inGeoFile = "Input Geometry File";
static constexpr auto outMeshFile = "Output Mesh File";
static constexpr auto bmParams = "blockMesh Parameters";
static constexpr auto smParams = "snappyHexMesh Parameters";

}  // namespace JSON
}  // namespace DRV
}  // namespace NEM

NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_BASE(
    NEM::DRV::NemDriver,
    (IF_EPIC(NEM::DRV::InputGenDriver)
     IF_HDF5_GMSH(NEM::DRV::ProteusDriver)
     IF_CGNS(NEM::DRV::RocPartCommGenDriver)
     IF_TEMPLATE_MESH(NEM::DRV::TemplateMeshDriver)
     IF_OCC(NEM::DRV::NucMeshDriver)
     NEM::DRV::AutoVerificationDriver, NEM::DRV::ConversionDriver,
     NEM::DRV::MeshGenDriver, NEM::DRV::MeshQualityDriver,
     NEM::DRV::PackMeshDriver, NEM::DRV::RefineDriver,
     NEM::DRV::TransferDriver),
    1,
    (getProgramType, , NEM::DRV::JSON::programType, JSONCONS_RDONLY,
     NEM_JSON_CHECK_KEY_ONLY))

NEM_JSON_WRAP_SMART_PTR(NEM::DRV::NemDriver, std::unique_ptr)
NEM_JSON_WRAP_SMART_PTR(NEM::DRV::NemDriver, std::shared_ptr)

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::DriverOutFile,
                                (outputFile, NEM::DRV::JSON::outMeshFile))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::DriverInOutFiles,
                                (inputMeshFile, NEM::DRV::JSON::inMeshFile),
                                (outputMeshFile, NEM::DRV::JSON::outMeshFile))

namespace jsoncons {
template <typename Json>
struct json_type_traits<Json, void *> {
  static constexpr bool is(const Json &ajson) { return false; }
  static void *as(const Json &ajson) { return nullptr; }
  static Json to_json(const void *&aval) { return Json::null(); }
};
template <>
struct is_json_type_traits_declared<void *> : public std::true_type {};
}  // namespace jsoncons

#undef IF_CGNS
#undef IF_EPIC
#undef IF_HDF5
#undef IF_TEMPLATE_MESH
#undef IF_OCC

#ifdef HAVE_EPIC
#  include "Drivers/InputGenJson.H"
#endif
#if defined(HAVE_HDF5) && defined(HAVE_GMSH)
#  include "Drivers/ProteusJson.H"
#endif
#ifdef HAVE_CGNS
#  include "Drivers/RocPartCommGenJson.H"
#endif
#ifdef HAVE_TEMPLATE_MESH
#  include "Drivers/TemplateMeshJson.H"
#endif
#ifdef HAVE_OCC
#  include "Drivers/NucMeshJson.H"
#endif
#include "Drivers/AutoVerificationJson.H"
#include "Drivers/Conversion/ConversionJson.H"
#include "Drivers/MeshGen/MeshGenJson.H"
#include "Drivers/MeshQualityJson.H"
#include "Drivers/PackMesh/PackMeshJson.H"
#include "Drivers/Refine/RefineJson.H"
#include "Drivers/TransferJson.H"

#endif  // NEMOSYS_DRIVERJSONTYPETRAITS_H_
