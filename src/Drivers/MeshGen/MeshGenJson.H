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
#ifndef NEMOSYS_MESHGENJSON_H_
#define NEMOSYS_MESHGENJSON_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <jsoncons/json.hpp>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/MeshGen/MeshGenDriver.H"

#ifdef HAVE_CFMSH
#  define HAVE_MSHGEN
#  define IF_CFMSH(a) a,
#  include "Drivers/MeshGen/BlockMeshMeshGenDriver.H"
#  include "Drivers/MeshGen/CFMeshMeshGenDriver.H"
#  include "Drivers/MeshGen/SnappyMeshMeshGenDriver.H"
#else
#  define IF_CFMSH(a)
#endif

#ifdef HAVE_NGEN
#  define HAVE_MSHGEN
#  define IF_NGEN(a) a,
#  include "Drivers/MeshGen/NetgenMeshGenDriver.H"
#else
#  define IF_NGEN(a)
#endif

#ifdef HAVE_GMSH
#  define HAVE_MSHGEN
#  define IF_GMSH(a) a,
#  include "Drivers/MeshGen/GmshMeshGenDriver.H"
#else
#  define IF_GMSH(a)
#endif

#ifndef HAVE_MSHGEN
#  error "No mesh generators enabled. Check configuration settings."
#endif

namespace NEM {
namespace DRV {
namespace JSON {

static constexpr auto meshGenOpts = "Mesh Generation Options";
static constexpr auto meshGenEngine = "Mesh Generation Engine";

struct colorMapPair {
  std::array<int, 3> color;
  std::string group;
};

inline std::vector<colorMapPair> gmshColorMapToJSON(
    const std::map<std::array<int, 3>, std::string> &colorMap) {
  std::vector<colorMapPair> out(colorMap.size());
  std::transform(colorMap.begin(), colorMap.end(), out.begin(),
                 [](const std::pair<const std::array<int, 3>, std::string> &x) {
                   return colorMapPair{x.first, x.second};
                 });
  return out;
}

inline std::map<std::array<int, 3>, std::string> gmshColorMapFromJSON(
    const std::vector<colorMapPair> &colorMap) {
  std::map<std::array<int, 3>, std::string> out;
  for (const auto &color : colorMap) {
    out.emplace(color.color, color.group);
  }
  return out;
}

}  // namespace JSON
}  // namespace DRV
}  // namespace NEM

// MeshGenDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_INTERMEDIATE(
    NEM::DRV::MeshGenDriver,
    (IF_CFMSH(NEM::DRV::BlockMeshMeshGenDriver)
     IF_CFMSH(NEM::DRV::CFMeshMeshGenDriver)
     IF_CFMSH(NEM::DRV::SnappyMeshMeshGenDriver)
     IF_NGEN(NEM::DRV::NetgenMeshGenDriver)
     IF_GMSH(NEM::DRV::GmshMeshGenDriver)
     void),
    NEM::DRV::NemDriver, 2,
    (, , NEM::DRV::JSON::meshFiles, NEM_JSON_RDONLY_OVERRIDE,
     NEM_JSON_CHECK_KEY_ONLY),
    (, , NEM::DRV::JSON::meshGenOpts, NEM_JSON_RDONLY_OVERRIDE,
     NEM_JSON_CHECK_KEY_ONLY),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::MeshGenDriver::programType;
     }))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::MeshGenDriver::MeshGenFiles,
                                (inputGeoFile, NEM::DRV::JSON::inGeoFile),
                                (outputMeshFile, NEM::DRV::JSON::outMeshFile))

#ifdef HAVE_CFMSH
// BlockMeshMeshGenDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::BlockMeshMeshGenDriver, NEM::DRV::MeshGenDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, NEM::DRV::JSON::meshGenOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::BlockMeshMeshGenDriver::Opts,
    (meshEngine, NEM::DRV::JSON::meshGenEngine, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::BlockMeshMeshGenDriver::Opts::meshEngine;
     }),
    (params, NEM::DRV::JSON::bmParams))

JSONCONS_N_MEMBER_NAME_TRAITS(
    blockMeshParams, 1, (cnvrtToMeters, "scaleToMeters"),
    (ownBlockMshDict, "Input Dict File"), (cellSize, "Cell_Size"),
    (nCells, "NumCells"),
    (shape, "Block Parameters", JSONCONS_RDWR, NEM_JSON_CHECK_KEY_ONLY,
     std::dynamic_pointer_cast<bmBox>,
     (std::dynamic_pointer_cast<bmShape, bmBox>)),
    (shape, "Sphere Parameters", JSONCONS_RDWR, NEM_JSON_CHECK_KEY_ONLY,
     std::dynamic_pointer_cast<bmSphere>,
     (std::dynamic_pointer_cast<bmShape, bmSphere>)),
    (shape, "Cylinder/Tapered_Cone Parameters", JSONCONS_RDWR,
     NEM_JSON_CHECK_KEY_ONLY, std::dynamic_pointer_cast<bmCylTaperedCone>,
     (std::dynamic_pointer_cast<bmShape, bmCylTaperedCone>)))

JSONCONS_N_MEMBER_NAME_TRAITS(bmBox, 0, (autoGenerate, "Auto_Generate"),
                              (init, "InitialPoint"), (len, "Length"),
                              (smplGrading, "Grading"))

NEM_JSON_SMART_PTR_VAL(bmBox, std::shared_ptr)

JSONCONS_N_MEMBER_NAME_TRAITS(bmAutoGenBox, 0, (packFileName, "InputPackFile"),
                              (offset, "Offset"))

JSONCONS_N_MEMBER_NAME_TRAITS(bmSphere, 2, (center, "Center"),
                              (radius, "Radius"), (sphrGrading, "Grading"))

NEM_JSON_SMART_PTR_VAL(bmSphere, std::shared_ptr)

JSONCONS_N_MEMBER_NAME_TRAITS(bmCylTaperedCone, 3, (centerCyl, "Center"),
                              (height, "Height"), (radius1, "Radius1"),
                              (radius2, "Radius2"), (cylGrading, "Grading"))

NEM_JSON_SMART_PTR_VAL(bmCylTaperedCone, std::shared_ptr)

// CFMeshMeshGenDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::CFMeshMeshGenDriver, NEM::DRV::MeshGenDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, NEM::DRV::JSON::meshGenOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::CFMeshMeshGenDriver::Opts,
    (meshEngine, NEM::DRV::JSON::meshGenEngine, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::CFMeshMeshGenDriver::Opts::meshEngine;
     }),
    (params, "CFMesh Parameters"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    cfmeshParams, 1, (generator, "Generator"), (maxCellSize, "MaxCellSize"),
    (minCellSize, "MinCellSize"), (bndryCellSize, "BoundaryCellSize"),
    (bndryCellSizeRefThk, "BoundaryCellSizeRefinementThickness"),
    (keepCellIB, "KeepCellsIntersectingBoundary"),
    (chkGluMsh, "CheckForGluedMesh"),
    (alwDiscDomains, "AllowDisconnectedDomains"),
    (boundaryLayers, "BoundaryLayers"), (srfEdge, "SurfaceFeatureEdges"),
    (objRefLst, "ObjectRefinements"),
    (improveMeshQuality, "ImproveMeshQuality"), (refPatches, "LocalRefinement"),
    (renBndry, "RenameBoundary"))

JSONCONS_N_MEMBER_NAME_TRAITS(cfmBoundaryLayer, 2, (blNLyr, "NLayers"),
                              (blThkRto, "ThicknessRatio"),
                              (maxFrstLyrThk, "MaxFirstLayerThickness"),
                              (alwDiscont, "AllowDiscontinuity"),
                              (blPatches, "PatchBoundaryLayers"))

JSONCONS_N_MEMBER_NAME_TRAITS(cfmPtchBndLyr, 1, (patchName, "PatchName"),
                              (alwDiscont, "AllowDiscontinuity"),
                              (maxFrstLyrThk, "MaxFirstLayerThickness"),
                              (blNLyr, "NLayers"), (blThkRto, "ThicknessRatio"))

JSONCONS_N_MEMBER_NAME_TRAITS(cfmSrfFeatEdge, 0, (srfEdgAng, "Angle"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(cfmObjRef, (name, "Name"), (params, "Params"))

JSONCONS_N_MEMBER_NAME_TRAITS(cfmMeshQual, 4, (qltNItr, "NIterations"),
                              (qltNLop, "NLoops"),
                              (qltQltThr, "QualityThreshold"),
                              (qltNSrfItr, "NSurfaceIterations"),
                              (qltConCelSet, "ConstrainedCellsSet"))

JSONCONS_N_MEMBER_NAME_TRAITS(cfmLclRefPatch, 1, (patchName, "PatchName"),
                              (aditRefLvls, "AdditionalRefinementLevels"),
                              (refThickness, "RefinementThickness"),
                              (cellSize, "CellSize"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(cfmRenBndry, (defName, "DefaultName"),
                                (defType, "DefaultType"),
                                (newPatches, "NewPatchNames"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(cfmNewPatch, (name, "Name"),
                                (newName, "NewName"), (newType, "NewType"))

// SnappyMeshMeshGenDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::SnappyMeshMeshGenDriver, NEM::DRV::MeshGenDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, NEM::DRV::JSON::meshGenOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::SnappyMeshMeshGenDriver::Opts,
    (meshEngine, NEM::DRV::JSON::meshGenEngine, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::SnappyMeshMeshGenDriver::Opts::meshEngine;
     }),
    (params, NEM::DRV::JSON::smParams))

JSONCONS_N_MEMBER_NAME_TRAITS(
    snappymeshParams, 3, (withCastMesh, "Castellated Mesh"),
    (withSnap, "Snapping"), (withLayers, "Layer Addition"),
    (mergeTol, "mergeTolerance"), (geoDef, "Geometry Definition"),
    (castMeshControls, "Castellated Mesh Controls"),
    (snapControls, "Snapping Controls"), (layerControls, "Layer Controls"),
    (qualityControls, "Mesh Quality Controls"))

JSONCONS_N_MEMBER_NAME_TRAITS(shmGeo, 1,
                              (withMultiPatches, "Enable Multi Patches"),
                              (singleSolidPatch, "InputPatchName"),
                              (stlPatchDefs, "Geometry Patches"),
                              (srchShape, "Custom Patches"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(shmSTLDefinition,
                                (STLPatchName, "Geometry Patch Name"),
                                (snappyPatchName, "Output Patch Name"))

NEM_JSON_N_MEMBER_NAME_TRAITS_BASE(shmSearchableShape,
                                   (shmSearchableBox, shmSearchableCylinder,
                                    shmSearchableSphere),
                                   2, (patchName, "Custom Patch Name"),
                                   (getType(), "Searchable Shape",
                                    JSONCONS_RDONLY))

NEM_JSON_WRAP_SMART_PTR(shmSearchableShape, std::shared_ptr)

NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(shmSearchableBox, shmSearchableShape, 2,
                                    (minBound, "minimum bound"),
                                    (maxBound, "maximum bound"))

NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(shmSearchableCylinder, shmSearchableShape,
                                    3, (axisPoint1, "Axis Point 1"),
                                    (axisPoint2, "Axis Point 2"),
                                    (radius, "Radius"))

NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(shmSearchableSphere, shmSearchableShape, 2,
                                    (center, "Center"), (radius, "Radius"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    shmCastMeshControls, 1, (withCellZones, "CellZones"),
    (castMeshGpLvl, "GeneralGapLevelIncrement"), (maxLCells, "maxLocalCells"),
    (maxGCells, "maxGlobalCells"), (minRefCells, "minRefCells"),
    (cellsBetnLvls, "nCellsBetweenLevels"),
    (refSurfLvlMin, "surfaceRefinementLvlMin"),
    (refSurfLvlMax, "surfaceRefinementLvlMax"),
    (featAngle, "resolveFeatureAngle"), (gapLvlInc, "gapLevelIncrement"),
    (planarAngle, "planarAngle"), (locMesh, "locationInMesh"),
    (alwFreeZone, "allowFreeStandingZoneFaces"), (ftrEdge, "Feature File"),
    (surfRefs, "SurfaceRefinementRegions"), (geomRefs, "GeomRefinementRegions"))

JSONCONS_N_MEMBER_NAME_TRAITS(shmFeatureEdgeRef, 2, (fileName, "File Name"),
                              (maxLvl, "MaxLevel"), (minLvl, "MinLevel"))

JSONCONS_N_MEMBER_NAME_TRAITS(shmSurfRefine, 2, (refPatchName, "Patch Name"),
                              (maxLvl, "MaxLevel"), (patchType, "Patch Type"),
                              (minLvl, "MinLevel"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(shmRegionRefine, (patchName, "Patch Name"),
                                (mode, "Mode"), (minLvl, "MinLevel"),
                                (maxLvl, "MaxLevel"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    shmSnapControls, 0, (nSmoothPatch, "nSmoothPatch"),
    (tolerance, "tolerance"), (nSolveIter, "snapSolveIter"),
    (nRelaxIter, "snapRelaxIter"), (nFeatureSnapIter, "nFeatureSnapIter"),
    (implicitFeatureSnap, "implicitFeatureSnap"),
    (explicitFeatureSnap, "explicitFeatureSnap"),
    (multiRegionFeatureSnap, "multiRegionFeatureSnap"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    shmLayerControls, 0, (relSize, "relativeSizes"),
    (expRatio, "expansionRatio"), (finLThick, "finalLayerThickness"),
    (minThick, "minThickness"), (firstLyrThickness, "firstLayerThickness"),
    (thickness, "thickness"), (nGrow, "nGrow"), (featAngle, "featureAngle"),
    (relaxIter, "nRelaxIter"), (smthSurfNorm, "nSmoothSurfaceNormals"),
    (smthNorm, "nSmoothNormals"), (smthThick, "nSmoothThickness"),
    (maxFcTR, "maxFaceThicknessRatio"),
    (maxThickTMR, "maxThicknessToMedialRatio"),
    (minMedAngl, "minMedialAxisAngle"), (bufferCells, "nBufferCellsNoExtrude"),
    (nIter, "nLayerIter"), (nRelaxedIter, "nRelaxedIter"),
    (slipFeatureAngle, "slipFeatureAngle"),
    (nMedialAxisIter, "nMedialAxisIter"),
    (nSmoothDisplacement, "nSmoothDisplacement"), (layerVec, "Layers"))

JSONCONS_N_MEMBER_NAME_TRAITS(shmLayers, 1, (patchName, "Patch Name"),
                              (nSurfaceLayers, "nSurfaceLayers"),
                              (expansionRatio, "expansionRatio"),
                              (finalLayerThickness, "finalLayerThickness"),
                              (firstLyrThickness, "firstLyrThickness"),
                              (thickness, "thickness"),
                              (minThickness, "minThickness"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    shmMeshQualityControls, 0, (maxNonOrtho, "maxNonOrtho"),
    (maxBndrySkew, "maxBoundarySkewness"), (maxIntSkew, "maxInternalSkewness"),
    (maxConc, "maxConcave"), (minVol, "minVol"), (minTetQ, "minTetQuality"),
    (minArea, "minArea"), (minTwist, "minTwist"), (minFaceW, "minFaceWeight"),
    (minVolRto, "minVolRatio"), (minDet, "minDeterminant"),
    (minTriTwist, "minTriangleTwist"), (smoothScale, "qcnSmoothScale"),
    (errReduction, "errorReduction"))
#endif

#ifdef HAVE_NGEN
// NetgenMeshGenDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::NetgenMeshGenDriver, NEM::DRV::MeshGenDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, NEM::DRV::JSON::meshGenOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::NetgenMeshGenDriver::Opts,
    (meshEngine, NEM::DRV::JSON::meshGenEngine, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::NetgenMeshGenDriver::Opts::meshEngine;
     }),
    (params, "Netgen Parameters"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    netgenParams, 0, (uselocalh, "uselocalh"), (maxh, "maxh"),
    (fineness, "fineness"), (grading, "grading"),
    (elementsperedge, "elementsperedge"),
    (elementspercurve, "elementspercurve"),
    (closeedgeenable, "closeedgeenable"), (closeedgefact, "closeedgefact"),
    (second_order, "second_order"), (meshsize_filename, "meshsize_filename"),
    (quad_dominated, "quad_dominated"), (optvolmeshenable, "optvolmeshenable"),
    (optsteps_2d, "optsteps_2d"), (optsteps_3d, "optsteps_3d"),
    (invert_tets, "invert_tets"), (invert_trigs, "invert_trigs"),
    (check_overlap, "check_overlap"),
    (check_overlapping_boundary, "check_overlapping_boundary"),
    (refine_with_geom, "refine_with_geometry_adaptation"),
    (refine_without_geom, "refine_without_geometry_adaptation"))
#endif

#ifdef HAVE_GMSH
// GmshMeshGenDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::GmshMeshGenDriver, NEM::DRV::MeshGenDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, NEM::DRV::JSON::meshGenOpts, JSONCONS_RDWR))

JSONCONS_ALL_MEMBER_NAME_TRAITS(
    NEM::DRV::GmshMeshGenDriver::Opts,
    (meshEngine, NEM::DRV::JSON::meshGenEngine, JSONCONS_RDONLY,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::GmshMeshGenDriver::Opts::meshEngine;
     }),
    (params, "Gmsh Parameters"))

JSONCONS_N_MEMBER_NAME_TRAITS(
    NEM::GEN::gmshParams, 0, (minSize, "minSize"), (maxSize, "maxSize"),
    (algo2D, "surfaceAlgorithm"), (algo3D, "volumeAlgorithm"),
    (extSizeFromBoundary, "extendSizeFromBoundary"),
    (sizeFromCurvature, "sizeFromCurvature"),
    (minElePer2Pi, "minElementsPer2Pi"), (optimize, "optimize"),
    (optimizeThreshold, "optimizeThreshold"), (sizeFields, "SizeFields"),
    (bgField, "BackgroundField"), (elementOrder, "elementOrder"),
    (subdivisionAlg, "subdivisionAlgorithm"), (saveAll, "saveAll"),
    (fragmentAll, "fragment"),
    (color2groupMap, "ColorMap", JSONCONS_RDWR, NEM_JSON_CHECK_KEY_ONLY,
     NEM::DRV::JSON::gmshColorMapToJSON, NEM::DRV::JSON::gmshColorMapFromJSON),
    (transfiniteBlocks, "TransfiniteBlocks"))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::JSON::colorMapPair, (color, "Color"),
                                (group, "Group"))

namespace jsoncons {
template <typename Json>
struct json_type_traits<Json, NEM::GEN::TransfiniteBlock> {
  using value_type = NEM::GEN::TransfiniteBlock;
  using allocator_type = typename Json::allocator_type;
  static constexpr auto idStr = "Volume";
  static constexpr auto axisStr = "Axis";
  static constexpr auto vertStr = "Vertices";
  static constexpr auto progStr = "Progression";
  static constexpr auto bumpStr = "Bump";
  static bool is(const Json &ajson) noexcept {
    return ajson.is_object() && ajson.contains(idStr) &&
           ajson.contains(axisStr) && ajson.contains("x") &&
           ajson.contains("y") && ajson.contains("z");
  }
  static value_type as(const Json &ajson) {
    static const std::array<std::pair<jsoncons::string_view, int>, 3> axisMap{
        {{"x", 0}, {"y", 1}, {"z", 2}}};
    if (!is(ajson))
      throw convert_error(convert_errc::conversion_failed,
                          "Not a NEM::GEN::TransfiniteBlock");
    value_type aval{};
    json_traits_helper<Json>::set_udt_member(ajson, idStr, aval.id);
    json_traits_helper<Json>::set_udt_member(ajson, axisStr, aval.axis);
    for (int i = 0; i < 3; ++i) {
      auto &axis = aval.axis[i];
      double axis_len =
          std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
      if (axis_len < 1e-3) {
        std::cerr
            << "Axis " << i << " for block " << aval.id
            << " is too small. Please prescribe an axis with length > 1e-3"
            << std::endl;
        exit(1);
      }
      for (auto &val : axis) {
        val /= axis_len;
      }
    }
    for (const auto &axis : axisMap) {
      const auto &axisJson = ajson[axis.first];
      aval.vert[axis.second] = axisJson[vertStr].template as<int>();
      if (axisJson.contains(bumpStr)) {
        aval.type[axis.second] = bumpStr;
        aval.coef[axis.second] = axisJson[bumpStr].as_double();
      } else if (axisJson.contains(progStr)) {
        aval.type[axis.second] = progStr;
        aval.coef[axis.second] = axisJson[progStr].as_double();
      }
    }
    return aval;
  }
  static Json to_json(const value_type &aval,
                      allocator_type alloc = allocator_type()) {
    static const std::array<std::pair<jsoncons::string_view, int>, 3> axisMap{
        {{"x", 0}, {"y", 1}, {"z", 2}}};
    Json ajson(json_object_arg, semantic_tag::none, alloc);
    ajson.try_emplace(idStr, aval.id);
    ajson.try_emplace(axisStr, aval.axis);
    for (const auto &axis : axisMap) {
      auto &axisJson = ajson[axis.first];
      axisJson.try_emplace(vertStr, aval.vert[axis.second]);
      axisJson.try_emplace(aval.type[axis.second], aval.coef[axis.second]);
    }
    return ajson;
  }
};
template <>
struct is_json_type_traits_declared<NEM::GEN::TransfiniteBlock>
    : public std::true_type {};
}  // namespace jsoncons

namespace jsoncons {
template <typename Json>
struct json_type_traits<Json, NEM::GEN::volSizeField> {
  using value_type = NEM::GEN::volSizeField;
  using allocator_type = typename Json::allocator_type;
  static bool is(const Json &ajson) noexcept {
    return ajson.is_object() && ajson.contains("Type") &&
           ajson.contains("ID") && ajson.contains("Params");
  }
  static value_type as(const Json &ajson) {
    if (!is(ajson))
      throw convert_error(convert_errc::conversion_failed,
                          "Not a NEM::GEN::volSizeField");
    value_type aval{};
    json_traits_helper<Json>::set_udt_member(ajson, "Type", aval.type);
    json_traits_helper<Json>::set_udt_member(ajson, "ID", aval.id);
    for (const auto &param : ajson.at("Params").object_range()) {
      const auto &val = param.value();
      if (val.is_array()) {
        if (val.at(0).is_string()) {
          aval.strg_list_params.emplace_back(
              param.key(), val.template as<std::vector<std::string>>());
        } else {
          aval.num_list_params.emplace_back(
              param.key(), val.template as<std::vector<double>>());
        }
      } else {
        aval.params.emplace_back(param.key(), val.as_double());
      }
    }
    return aval;
  }
  static Json to_json(const value_type &aval,
                      allocator_type alloc = allocator_type()) {
    Json ajson(json_object_arg, semantic_tag::none, alloc);
    ajson.try_emplace("Type", aval.type);
    ajson.try_emplace("ID", aval.id);
    auto &params = ajson["Params"];
    for (const auto &param : aval.strg_list_params) {
      params.try_emplace(param.first, param.second);
    }
    for (const auto &param : aval.num_list_params) {
      params.try_emplace(param.first, param.second);
    }
    for (const auto &param : aval.params) {
      params.try_emplace(param.first, param.second);
    }
    return ajson;
  }
};
template <>
struct is_json_type_traits_declared<NEM::GEN::volSizeField>
    : public std::true_type {};
}  // namespace jsoncons
#endif

#undef IF_CFMSH
#undef IF_NGEN
#undef IF_GMSH
#undef HAVE_MSHGEN

#endif  // NEMOSYS_MESHGENJSON_H_
