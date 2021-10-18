#ifndef NEMOSYS_PACKMESHJSON_H_
#define NEMOSYS_PACKMESHJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/PackMesh/PackMeshDriver.H"

#ifdef HAVE_CFMSH
#  define IF_CFMSH(a) a,
#  include "Drivers/PackMesh/HexPackMeshDriver.H"
#else
#  define IF_CFMSH(a)
#endif

#ifdef HAVE_GMSH
#  define IF_GMSH(a) a,
#  include "Drivers/PackMesh/SurfacePackMeshDriver.H"
#else
#  define IF_GMSH(a)
#endif

namespace NEM {
namespace DRV {
namespace JSON {

static constexpr auto packMeshOpts = "Pack Mesh Options";

}  // namespace JSON
}  // namespace DRV
}  // namespace NEM

// PackMeshDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_INTERMEDIATE(
    NEM::DRV::PackMeshDriver,
    (IF_CFMSH(NEM::DRV::HexPackMeshDriver)
     IF_GMSH(NEM::DRV::SurfacePackMeshDriver)
     void),
    NEM::DRV::NemDriver, 1,
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::PackMeshDriver::programType;
     }))

#ifdef HAVE_GMSH
// SurfacePackMeshDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::SurfacePackMeshDriver, NEM::DRV::PackMeshDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, NEM::DRV::JSON::packMeshOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::SurfacePackMeshDriver::Files,
                                (rocpackFile, "Input Rocpack File"),
                                (outputMeshFile, NEM::DRV::JSON::outMeshFile))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::SurfacePackMeshDriver::Opts, 3,
    (type, "Type", JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::SurfacePackMeshDriver::Opts::type;
     }),
    (engine, "Engine", JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::SurfacePackMeshDriver::Opts::engine;
     }),
    (periodic3DOpts, "Periodic 3D Mesh"), (meshSize, "Mesh Size"),
    (meshAlgorithm, "Mesh Algorithm"), (scaleValue, "Scale Value"),
    (removeBoundaryPacks, "Remove geometries on boundary"),
    (enableDefaultOut, "Enalbe Default Outputs"),
    (preserveSize, "Enable Size Preservation"),
    (refineLevel, "Refinement Levels"), (upperThreshold, "Upper Threshold"),
    (lowerThreshold, "Lower Threshold"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::SurfacePackMeshDriver::Periodic3DOpts,
                              0, (physGrpOptions, "Physical Group Options"),
                              (createCohesive, "Create cohesive elements"),
                              (enablePatches, "Enable Patches"),
                              (setPeriodicGeo, "Set Periodic Geometry"),
                              (elemOrder, "Element Order"),
                              (customDomain, "Custom Domain"),
                              (transferMesh, "TransferMesh"))

JSONCONS_ENUM_NAME_TRAITS(NEM::DRV::SurfacePackMeshDriver::PhysGrpOpts,
                          (NONE, "None"), (MULTI, "Multi Physical Groups"),
                          (TWO, "Two Physical Groups"),
                          (PER_SHAPE, "Physical Group per Shape"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::SurfacePackMeshDriver::CustomDomain,
                                (initial, "Initial"), (length, "Length"))
#endif

#ifdef HAVE_CFMSH
// HexPackMeshDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::HexPackMeshDriver, NEM::DRV::PackMeshDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, NEM::DRV::JSON::packMeshOpts, JSONCONS_RDWR))

namespace jsoncons {
template <typename Json>
struct json_type_traits<Json, NEM::DRV::HexPackMeshDriver::Files> {
  using value_type = NEM::DRV::HexPackMeshDriver::Files;
  using allocator_type = typename Json::allocator_type;
  static constexpr auto rockFile = "Input Rocpack File";
  static constexpr auto pmFile = "Output Pack Mesh File";
  static constexpr auto smFile = "Output Surrounding Mesh File";
  static constexpr auto cmFile = "Output Combined Mesh File";
  static bool is(const Json &ajson) noexcept {
    return ajson.is_object() &&
           (ajson.contains(rockFile) !=
            ajson.contains(NEM::DRV::JSON::inGeoFile)) &&
           ajson.contains(pmFile) && ajson.contains(smFile);
  }
  static value_type as(const Json &ajson) {
    if (!is(ajson))
      throw convert_error(convert_errc::conversion_failed,
                          "Not a HexPackMeshDriver::Files");
    value_type aval{};
    if (ajson.contains(rockFile))
      aval.setRocpackFile(ajson.at(rockFile).as_string());
    if (ajson.contains(NEM::DRV::JSON::inGeoFile))
      aval.setGeoFile(ajson.at(NEM::DRV::JSON::inGeoFile).as_string());
    json_traits_helper<Json>::set_udt_member(ajson, pmFile,
                                             aval.outPackMeshFile);
    json_traits_helper<Json>::set_udt_member(ajson, smFile,
                                             aval.outSurroundingFile);
    if (ajson.contains(cmFile))
      json_traits_helper<Json>::set_udt_member(ajson, cmFile,
                                               aval.outCombinedFile);
    return aval;
  }
  static Json to_json(const value_type &aval,
                      allocator_type alloc = allocator_type()) {
    Json ajson(json_object_arg, semantic_tag::none, alloc);
    ajson.try_emplace(
        aval.isInputRocpackFile() ? rockFile : NEM::DRV::JSON::inGeoFile,
        aval.rocpackOrGeoFile);
    ajson.try_emplace(pmFile, aval.outPackMeshFile);
    ajson.try_emplace(smFile, aval.outSurroundingFile);
    ajson.try_emplace(cmFile, aval.outCombinedFile);
    return ajson;
  }
};
template <>
struct is_json_type_traits_declared<NEM::DRV::HexPackMeshDriver::Files>
    : public std::true_type {};
}  // namespace jsoncons

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::DRV::HexPackMeshDriver::Opts, 2,
    (type, "Type", JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::HexPackMeshDriver::Opts::type;
     }),
    (engine, "Engine", JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::HexPackMeshDriver::Opts::engine;
     }),
    (smParams, NEM::DRV::JSON::smParams), (bmParams, NEM::DRV::JSON::bmParams),
    (mmfMergeParams, "MeshManipulation mergeMeshes Parameters"),
    (mmfCreatePatchParams, "MeshManipulation createPatch Parameters"),
    (locAdjust, "locInMesh_adjust"))

JSONCONS_N_MEMBER_NAME_TRAITS(MeshManipulationFoamParams::MergeMeshes, 0,
                              (overwriteMergeMsh, "overwrite?"),
                              (masterCasePath, "Master Region Path"),
                              (addCasePath, "Add Region Path"))

JSONCONS_N_MEMBER_NAME_TRAITS(MeshManipulationFoamParams::CreatePatch, 0,
                              (surroundingName, "Surrounding PatchName"),
                              (packsName, "Packs PatchName"),
                              (srrndngPatchType, "Surrounding PatchType"),
                              (packsPatchType, "Packs PatchType"),
                              (overwritecpMsh, "overwrite?"))
#endif

#undef IF_CFMSH
#undef IF_GMSH

#endif  // NEMOSYS_PACKMESHJSON_H_
