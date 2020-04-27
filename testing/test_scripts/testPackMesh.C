#include <PackMeshDriver.H>
#include <gtest.h>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iterator>
#include <string>
#include "MeshManipulationFoam.H"
#include "MeshManipulationFoamParams.H"
#include "NemDriver.H"
#include "blockMeshGen.H"
#include "blockMeshParams.H"
#include "cfmeshGen.H"
#include "cfmeshParams.H"
#include "snappymeshGen.H"
#include "snappymeshParams.H"
#include "vtkMesh.H"

const char *inp_json;
meshBase *mesh;
meshBase *ref;
jsoncons::json inputjson;

// Aux functions
bool compareFiles(const std::string &p1, const std::string &p2) {
  std::ifstream f1(p1, std::ifstream::binary | std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary | std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false;  // file problem
  }

  if (f1.tellg() != f2.tellg()) {
    return false;  // size mismatch
  }

  // seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

// Test implementations
int generate(const char *jsonF) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson_tmp;
  inputStream >> inputjson_tmp;
  inputjson = inputjson_tmp[0];

  // TODO: refactor the input processing, use actual driver objtect
  //       intead of parsing here.
  std::string ifname =
      inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  std::string ofname1 =
      inputjson["Mesh File Options"]["Output Pack Mesh File"].as<std::string>();
  std::string ofname2 =
      inputjson["Mesh File Options"]["Output Surrounding Mesh File"]
          .as<std::string>();

  cfmeshParams *cfparams = new cfmeshParams();
  std::string defaults1 =
      inputjson["Meshing Parameters"]["CFMesh Parameters"].as<std::string>();
  jsoncons::json cfmparams =
      inputjson["Meshing Parameters"]["CFMesh Parameters"];

  // required params here
  // cad file
  if (inputjson["Mesh File Options"].contains("Input Geometry File"))
    cfparams->geomFilePath =
        inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else {
    std::cerr << "A geometry file should be supplied.\n";
    throw;
  }

  // mesh generator
  if (cfmparams.contains("Generator"))
    cfparams->generator = cfmparams["Generator"].as<std::string>();
  else {
    std::cerr << "A mesh generation method should be selected.\n";
    std::cerr << "Options: cartesian2D tetMesh\n";
    throw;
  }

  // rest of params are optional
  if (cfmparams.contains("MaxCellSize"))
    cfparams->maxCellSize = cfmparams["MaxCellSize"].as<double>();
  if (cfmparams.contains("MinCellSize"))
    cfparams->minCellSize = cfmparams["MinCellSize"].as<double>();
  if (cfmparams.contains("BoundaryCellSize"))
    cfparams->bndryCellSize = cfmparams["BoundaryCellSize"].as<double>();
  if (cfmparams.contains("KeepCellsIntersectingBoundary"))
    cfparams->keepCellIB =
        cfmparams["KeepCellsIntersectingBoundary"].as<double>();
  if (cfmparams.contains("CheckForGluedMesh"))
    cfparams->chkGluMsh = cfmparams["CheckForGluedMesh"].as<double>();
  if (cfmparams.contains("AllowDisconnectedDomains"))
    cfparams->_alwDiscDomains =
        cfmparams["AllowDisconnectedDomains"].as<bool>();

  // optional capability
  std::string cap = "BoundaryLayers";
  if (cfmparams.contains(cap)) {
    cfparams->_withBndLyr = true;
    cfparams->blNLyr = cfmparams[cap]["NLayers"].as<double>();
    cfparams->blThkRto = cfmparams[cap]["ThicknessRatio"].as<double>();
    if (cfmparams[cap].contains("MaxFirstLayerThickness"))
      cfparams->maxFrstLyrThk =
          cfmparams[cap]["MaxFirstLayerThickness"].as<double>();
    if (cfmparams[cap].contains("AllowDiscontinuity"))
      cfparams->alwDiscont = cfmparams[cap]["AllowDiscontinuity"].as<bool>();

    // patch boundary layers
    std::string subcap = "PatchBoundaryLayers";
    if (cfmparams[cap].contains(subcap)) {
      cfparams->_withBndLyrPtch = true;
      for (auto jptch : cfmparams[cap][subcap].array_range()) {
        cfmPtchBndLyr blPatch;
        blPatch.patchName = jptch["PatchName"].as<std::string>();
        if (jptch.contains("AllowDiscontinuity"))
          blPatch.alwDiscont = jptch["AllowDiscontinuity"].as<bool>();
        else
          blPatch.alwDiscont = false;
        if (jptch.contains("MaxFirstLayerThickness"))
          blPatch.maxFrstLyrThk = jptch["MaxFirstLayerThickness"].as<int>();
        else
          blPatch.maxFrstLyrThk = -1;
        if (jptch.contains("NLayers"))
          blPatch.blNLyr = jptch["NLayers"].as<int>();
        else
          blPatch.blNLyr = -1;
        if (jptch.contains("ThicknessRatio"))
          blPatch.blThkRto = jptch["ThicknessRatio"].as<double>();
        else
          blPatch.blThkRto = -1.;
        (cfparams->blPatches).push_back(blPatch);
      }
    }
  }

  // optional capability
  cap = "SurfaceFeatureEdges";
  if (cfmparams.contains(cap)) {
    cfparams->_withSrfEdg = true;
    cfparams->srfEdgAng = cfmparams[cap]["Angle"].as<double>();
  }

  // optional capability
  cap = "ObjectRefinements";
  if (cfmparams.contains(cap)) {
    cfparams->_withObjRfn = true;
    for (auto refObj : cfmparams[cap].array_range()) {
      cfmObjRef objRef;
      objRef.name = refObj["Name"].as<std::string>();
      for (const auto &prm : refObj["Params"].object_range()) {
        std::string key = std::string(prm.key());
        std::string val = prm.value().as<std::string>();
        objRef.params[key] = val;
      }
      (cfparams->objRefLst).push_back(objRef);
    }
  }

  // optional capability
  cap = "ImproveMeshQuality";
  if (cfmparams.contains(cap)) {
    cfparams->_withMshQlt = true;
    cfparams->qltNItr = cfmparams[cap]["NIterations"].as<int>();
    cfparams->qltNLop = cfmparams[cap]["NLoops"].as<int>();
    cfparams->qltQltThr = cfmparams[cap]["QualityThreshold"].as<double>();
    cfparams->qltNSrfItr = cfmparams[cap]["NSurfaceIterations"].as<int>();
    cfparams->qltConCelSet =
        cfmparams[cap].get_with_default("ConstrainedCellsSet", "none");
  }

  // optional capability
  cap = "LocalRefinement";
  if (cfmparams.contains(cap)) {
    cfparams->_withLclRef = true;
    for (auto jptch : cfmparams[cap].array_range()) {
      cfmLclRefPatch refPatch;
      refPatch.patchName = jptch["PatchName"].as<std::string>();
      if (jptch.contains("AdditionalRefinementLevels"))
        refPatch.aditRefLvls = jptch["AdditionalRefinementLevels"].as<int>();
      else
        refPatch.aditRefLvls = -1;
      if (jptch.contains("CellSize"))
        refPatch.cellSize = jptch["CellSize"].as<double>();
      else
        refPatch.cellSize = -1.;
      (cfparams->refPatches).push_back(refPatch);
    }
  }

  // optional capability
  cap = "RenameBoundary";
  if (cfmparams.contains(cap)) {
    cfparams->_withRenBndry = true;
    cfmRenBndry renBndry;
    renBndry.defName = cfmparams[cap]["DefaultName"].as<std::string>();
    renBndry.defType = cfmparams[cap]["DefaultType"].as<std::string>();
    for (auto jnw : cfmparams[cap]["NewPatchNames"].array_range()) {
      cfmNewPatch nwPatch = std::make_tuple(jnw["Name"].as<std::string>(),
                                            jnw["NewName"].as<std::string>(),
                                            jnw["NewType"].as<std::string>());
      renBndry.newPatches.push_back(nwPatch);
    }
    (cfparams->renBndry) = renBndry;
  }

  // snappy
  snappymeshParams *params = new snappymeshParams();
  std::string defaults2 =
      inputjson["Meshing Parameters"]["snappyHexMesh Parameters"]
          .as<std::string>();
  jsoncons::json shmparams =
      inputjson["Meshing Parameters"]["snappyHexMesh Parameters"];

  if (inputjson["Mesh File Options"].contains("Input Geometry File"))
    params->geomFileName =
        inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else {
    std::cerr << "A geometry file should be supplied.\n";
    throw;
  }

  jsoncons::json geomParams = shmparams["Geometry Definition"];
  jsoncons::json castMeshParams = shmparams["Castellated Mesh Controls"];
  jsoncons::json snapParams = shmparams["Snapping Controls"];
  jsoncons::json layerParams = shmparams["Mesh Layers Controls"];
  jsoncons::json qcMeshParams = shmparams["Mesh Quality Controls"];

  if (inputjson["Mesh File Options"].contains("Input Geometry File"))
    params->geomFileName =
        inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else {
    params->geomFileName = ifname + ".stl";
  }

  // General booleans
  if (shmparams.contains("Castellated Mesh"))
    params->_withCastMesh = shmparams["Castellated Mesh"].as<bool>();
  else {
    std::cerr << "Please specify on/off choice using \"Castellated Mesh\""
              << "\n"
              << std::endl;
    throw;
  }
  if (shmparams.contains("Snapping"))
    params->_withSnap = shmparams["Snapping"].as<bool>();
  else {
    std::cerr << "Please specify on/off choice using \"Snapping\""
              << "Keyword!\n"
              << std::endl;
    throw;
  }
  if (shmparams.contains("Layer Addition"))
    params->_withLayers = shmparams["Layer Addition"].as<bool>();
  else {
    std::cerr << "Please specify on/off choice using \"Layer Addition\""
              << "Keyword!\n"
              << std::endl;
    throw;
  }

  // Geometry definition with patch defining ability
  bool multiPatches = false;
  if (geomParams.contains("Enable Multi Patches"))
    multiPatches = geomParams["Enable Multi Patches"].as<bool>();
  else {
    std::cerr << "Please provide boolean \"Enable Multi Patches\".\n";
    throw;
  }

  if (!multiPatches) {
    if (geomParams.contains("InputPatchName"))
      params->singleSolidPatch = geomParams["InputPatchName"].as<std::string>();
    else {
      std::cerr << "A patch name for input geometry must be defined.\n";
      throw;
    }
  } else {
    // Geometry STL Patch Definition
    std::string cap4 = "Geometry Patches";
    if (geomParams.contains(cap4)) {
      params->_withMultiPatches = true;

      for (auto jptch3 : geomParams[cap4].array_range()) {
        shmSTLDefinition stlPatches;

        if (jptch3.contains("Geometry Patch Name"))
          stlPatches.STLPatchName =
              jptch3["Geometry Patch Name"].as<std::string>();
        else {
          std::cerr << "Please provide Output Patch Name for STL file provided"
                    << std::endl;
          throw;
        }

        if (jptch3.contains("Output Patch Name"))
          stlPatches.snappyPatchName =
              jptch3["Output Patch Name"].as<std::string>();
        else {
          std::cerr << "Please provide Output Patch Name for STL file provided"
                    << std::endl;
          throw;
        }

        (params->stlPatchDefs).push_back(stlPatches);
      }
    } else {
      std::cerr << "User has selected to define multiple geometry patches.\n"
                << "Please define them using \"Geometry Patches\" keyword"
                << std::endl;
      throw;
    }
  }

  // Define searchable shapes with whatever patch name user wants
  std::string cap5 = "Custom Patches";
  if (geomParams.contains(cap5)) {
    for (auto jptch3 : geomParams[cap5].array_range()) {
      shmSearchableShapes SSS;
      if (jptch3.contains("Custom Patch Name"))
        SSS.patchNm = jptch3["Custom Patch Name"].as<std::string>();
      else {
        std::cerr << "Please provide Custom Patch Name" << std::endl;
        throw;
      }

      if (jptch3.contains("Searchable Shape"))
        SSS.searchableName = jptch3["Searchable Shape"].as<std::string>();
      else {
        std::cerr << "Please provide Searchable Shape Name" << std::endl;
        throw;
      }

      if (SSS.searchableName == "searchableBox") {
        if (jptch3.contains("minimum bound"))
          SSS.shapeParameters1 = jptch3["minimum bound"].as<std::string>();
        else {
          std::cerr << "Please provide minimum bound (i.e (-1 -1 -1))"
                    << std::endl;
          throw;
        }
        if (jptch3.contains("maximum bound"))
          SSS.shapeParameters2 = jptch3["maximum bound"].as<std::string>();
        else {
          std::cerr << "Please provide maximum bound (i.e (1 1 1))"
                    << std::endl;
          throw;
        }
      } else if (SSS.searchableName == "searchableCylinder") {
        if (jptch3.contains("Axis Point 1"))
          SSS.shapeParameters1 = jptch3["Axis Point 1"].as<std::string>();
        else {
          std::cerr << "Please provide Axis Point 1 (i.e (-1 -1 -1))"
                    << std::endl;
          throw;
        }
        if (jptch3.contains("Axis Point 2"))
          SSS.shapeParameters2 = jptch3["Axis Point 2"].as<std::string>();
        else {
          std::cerr << "Please provide Axis Point 2 (i.e (1 1 1))" << std::endl;
          throw;
        }
        if (jptch3.contains("Radius"))
          SSS.rad = jptch3["Radius"].as<double>();
        else {
          std::cerr << "Please provide Radius" << std::endl;
          throw;
        }
      } else if (SSS.searchableName == "searchableSphere") {
        if (jptch3.contains("Center"))
          SSS.shapeParameters1 = jptch3["Center"].as<std::string>();
        else {
          std::cerr << "Please provide Center (i.e (1 1 1))" << std::endl;
          throw;
        }
        if (jptch3.contains("Radius"))
          SSS.rad = jptch3["Radius"].as<double>();
        else {
          std::cerr << "Please provide Radius" << std::endl;
          throw;
        }
      } else {
        std::cerr << SSS.searchableName << " is not supported yet!"
                  << std::endl;
        throw;
      }
      (params->srchShape).push_back(SSS);
    }
  }

  // Castellated Mesh Controls Inputs
  if (castMeshParams.contains("CellZones"))
    params->_withCellZones = castMeshParams["CellZones"].as<bool>();
  else {
    std::cerr << "Please specify on/off choice using \"CellZones\""
              << "Keyword!\n"
              << std::endl;
    throw;
  }
  if (castMeshParams.contains("RegionRefine"))
    params->_withGeomRefReg = castMeshParams["RegionRefine"].as<bool>();
  else {
    std::cerr << "Please specify on/off choice using \"RegionRefine\""
              << "Keyword!\n"
              << std::endl;
    throw;
  }
  if (castMeshParams.contains("SurfaceRefine"))
    params->_withSurfRefReg = castMeshParams["SurfaceRefine"].as<bool>();
  else {
    std::cerr << "Please specify on/off choice using \"SurfaceRefine\""
              << "Keyword!\n"
              << std::endl;
    throw;
  }
  if (castMeshParams.contains("GeneralGapLevelIncrement"))
    params->castMeshGpLvl =
        castMeshParams["GeneralGapLevelIncrement"].as<int>();
  else {
    params->castMeshGpLvl = 1;
  }
  if (castMeshParams.contains("maxLocalCells"))
    params->maxLCells = castMeshParams["maxLocalCells"].as<int>();
  else {
    params->maxLCells = 2000000;
  }
  if (castMeshParams.contains("maxGlobalCells"))
    params->maxGCells = castMeshParams["maxGlobalCells"].as<int>();
  else {
    params->maxGCells = 4000000;
  }
  if (castMeshParams.contains("minRefCells"))
    params->minRefCells = castMeshParams["minRefCells"].as<int>();
  else {
    params->minRefCells = 0;
  }
  if (castMeshParams.contains("nCellsBetweenLevels"))
    params->cellsBetnLvls = castMeshParams["nCellsBetweenLevels"].as<int>();
  else {
    params->cellsBetnLvls = 3;
  }
  if (castMeshParams.contains("surfaceRefinementLvlMin"))
    params->refSurfLvlMin = castMeshParams["surfaceRefinementLvlMin"].as<int>();
  else {
    params->refSurfLvlMin = 0;
  }
  if (castMeshParams.contains("surfaceRefinementLvlMax"))
    params->refSurfLvlMax = castMeshParams["surfaceRefinementLvlMax"].as<int>();
  else {
    params->refSurfLvlMax = 0;
  }
  if (castMeshParams.contains("resolveFeatureAngle"))
    params->featAngle = castMeshParams["resolveFeatureAngle"].as<double>();
  else {
    params->featAngle = 60;
  }
  if (castMeshParams.contains("gapLevelIncrement"))
    params->gPLvlInc = castMeshParams["gapLevelIncrement"].as<int>();
  else
    params->gPLvlInc = 1;
  if (castMeshParams.contains("planarAngle"))
    params->planarAngle = castMeshParams["planarAngle"].as<int>();
  else
    params->planarAngle = 1;
  if (castMeshParams.contains("locationInMeshX"))
    params->locMeshX = castMeshParams["locationInMeshX"].as<double>();
  else {
    std::cerr << "Location of a point in region you want to"
              << "keep cells is needed!" << std::endl;
    throw;
  }
  if (castMeshParams.contains("locationInMeshY"))
    params->locMeshY = castMeshParams["locationInMeshY"].as<double>();
  else {
    std::cerr << "Location of a point in region you want to"
              << "keep cells is needed!" << std::endl;
    throw;
  }
  if (castMeshParams.contains("locationInMeshZ"))
    params->locMeshZ = castMeshParams["locationInMeshZ"].as<double>();
  else {
    std::cerr << "Location of a point in region you want to"
              << "keep cells is needed!" << std::endl;
    throw;
  }
  if (castMeshParams.contains("allowFreeStandingZoneFaces"))
    params->_alwFreeZone =
        castMeshParams["allowFreeStandingZoneFaces"].as<bool>();
  else {
    params->_alwFreeZone = true;
  }

  // Refinement Surfaces
  std::string cap3 = "SurfaceRefinementRegions";
  if (castMeshParams.contains(cap3)) {
    params->_withSurfRefReg = true;

    for (auto jptch3 : castMeshParams[cap3].array_range()) {
      shmSurfRefine surfRef;

      if (jptch3.contains("Patch Name"))
        surfRef.refPatchNm = jptch3["Patch Name"].as<std::string>();
      else {
        surfRef.refPatchNm = params->singleSolidPatch;
      }

      if (jptch3.contains("Patch Type"))
        surfRef.patchType = jptch3["Patch Type"].as<int>();
      else {
        surfRef.patchType = "NO";
      }

      if (jptch3.contains("MinLevel"))
        surfRef.minLvl = jptch3["MinLevel"].as<int>();
      else {
        surfRef.minLvl = 1;
      }

      if (jptch3.contains("MaxLevel"))
        surfRef.maxLvl = jptch3["MaxLevel"].as<int>();
      else {
        std::cerr << "Please define maximum surface refinement"
                  << "\"MaxLevel\"\n"
                  << std::endl;
        throw;
      }
      (params->surfRefs).push_back(surfRef);
    }
  }

  // Features
  std::string cap6 = "Feature File";
  if (castMeshParams.contains(cap6)) {
    params->_withFeatureEdgeFile = true;
    for (auto jptch3 : castMeshParams[cap6].array_range()) {
      shmFeatureEdgeRef ftrOne;
      if (jptch3.contains("File Name"))
        ftrOne.fileName = jptch3["File Name"].as<int>();
      else {
        std::cerr << "Please provide feature file name" << std::endl;
        throw;
      }

      if (jptch3.contains("MinLevel"))
        ftrOne.minLvl = jptch3["MinLevel"].as<int>();
      else {
        ftrOne.minLvl = 1;
      }

      if (jptch3.contains("MaxLevel"))
        ftrOne.maxLvl = jptch3["MaxLevel"].as<int>();
      else {
        std::cerr << "Please define maximum surface refinement"
                  << "\"MaxLevel\"\n"
                  << std::endl;
        throw;
      }
      (params->ftrEdge).push_back(ftrOne);
    }
  }

  // Snapping Controls
  if (snapParams.contains("nSmoothPatch"))
    params->snapSmthPatch = snapParams["nSmoothPatch"].as<int>();
  else {
    params->snapSmthPatch = 4;
  }
  if (snapParams.contains("tolerance"))
    params->snapTol = snapParams["tolerance"].as<double>();
  else {
    params->snapTol = 0.5;
  }
  if (snapParams.contains("snapSolveIter"))
    params->solveSnapIter = snapParams["snapSolveIter"].as<int>();
  else {
    params->solveSnapIter = 200;
  }
  if (snapParams.contains("snapRelaxIter"))
    params->relaxSnapIter = snapParams["snapRelaxIter"].as<int>();
  else {
    params->relaxSnapIter = 6;
  }
  if (snapParams.contains("nFeatureSnapIter"))
    params->nFeatureSnapIter = snapParams["nFeatureSnapIter"].as<int>();
  else
    params->nFeatureSnapIter = 10;
  if (snapParams.contains("implicitFeatureSnap"))
    params->implicitFeatureSnap = snapParams["implicitFeatureSnap"].as<bool>();
  else
    params->implicitFeatureSnap = false;
  if (snapParams.contains("explicitFeatureSnap"))
    params->explicitFeatureSnap = snapParams["explicitFeatureSnap"].as<bool>();
  else
    params->explicitFeatureSnap = true;
  if (snapParams.contains("multiRegionFeatureSnap"))
    params->multiRegionFeatureSnap =
        snapParams["multiRegionFeatureSnap"].as<bool>();
  else
    params->multiRegionFeatureSnap = false;

  // Layer controls
  if (layerParams.contains("relativeSizes"))
    params->_relSize = layerParams["relativeSizes"].as<bool>();
  else {
    params->_relSize = 1;
  }
  if (layerParams.contains("expansionRatio"))
    params->expRatio = layerParams["expansionRatio"].as<double>();
  else {
    params->expRatio = 1.3;
  }
  if (layerParams.contains("finalLayerThickness"))
    params->finLThick = layerParams["finalLayerThickness"].as<double>();
  else {
    params->finLThick = 1.0;
  }
  if (layerParams.contains("minThickness"))
    params->minThick = layerParams["minThickness"].as<double>();
  else {
    params->minThick = 0.1;
  }
  if (layerParams.contains("nGrow"))
    params->nGrow = layerParams["nGrow"].as<int>();
  else {
    params->nGrow = 0;
  }
  if (layerParams.contains("featureAngle"))
    params->lyrFeatAngle = layerParams["featureAngle"].as<double>();
  else {
    params->lyrFeatAngle = 30;
  }
  if (layerParams.contains("nRelaxIter"))
    params->lyrRelaxIter = layerParams["nRelaxIter"].as<int>();
  else {
    params->lyrRelaxIter = 3;
  }
  if (layerParams.contains("nSmoothSurfaceNormals"))
    params->lyrSmthSurfNorm = layerParams["nSmoothSurfaceNormals"].as<int>();
  else {
    params->lyrSmthSurfNorm = 1;
  }
  if (layerParams.contains("nSmoothNormals"))
    params->lyrSmthNorm = layerParams["nSmoothNormals"].as<int>();
  else {
    params->lyrSmthNorm = 3;
  }
  if (layerParams.contains("nSmoothThickness"))
    params->lyrSmthThick = layerParams["nSmoothThickness"].as<int>();
  else {
    params->lyrSmthThick = 2;
  }
  if (layerParams.contains("maxFaceThicknessRatio"))
    params->lyrMaxFcTR = layerParams["maxFaceThicknessRatio"].as<double>();
  else {
    params->lyrMaxFcTR = 0.5;
  }
  if (layerParams.contains("maxThicknessToMedialRatio"))
    params->lyrMaxThickTMR =
        layerParams["maxThicknessToMedialRatio"].as<double>();
  else {
    params->lyrMaxThickTMR = 1.0;
  }
  if (layerParams.contains("minMedialAxisAngle"))
    params->lyrMinMedAngl = layerParams["minMedialAxisAngle"].as<double>();
  else {
    params->lyrMinMedAngl = 90;
  }
  if (layerParams.contains("nBufferCellsNoExtrude"))
    params->lyrBuffrCells = layerParams["nBufferCellsNoExtrude"].as<int>();
  else {
    params->lyrBuffrCells = 0;
  }
  if (layerParams.contains("nLayerIter"))
    params->lyrIter = layerParams["nLayerIter"].as<int>();
  else {
    params->lyrIter = 50;
  }
  if (layerParams.contains("nRelaxedIter"))
    params->nRelaxedIter = layerParams["nRelaxedIter"].as<int>();
  else {
    params->nRelaxedIter = 20;
  }
  if (layerParams.contains("slipFeatureAngle"))
    params->slipFeatureAngle = layerParams["slipFeatureAngle"].as<int>();
  else {
    params->slipFeatureAngle = 20;
  }

  // Layers
  std::string cap7 = "Layers";
  if (layerParams.contains(cap7)) {
    for (auto jptch3 : layerParams[cap7].array_range()) {
      shmLayers lyrOne;
      if (jptch3.contains("Patch Name"))
        lyrOne.patchName = jptch3["File Name"].as<int>();
      else {
        std::cerr << "Please provide patch name for layers" << std::endl;
        throw;
      }

      if (jptch3.contains("nSurfaceLayers"))
        lyrOne.nSurfaceLayers = jptch3["nSurfaceLayers"].as<int>();
      else {
        lyrOne.nSurfaceLayers = 1;
      }

      if (jptch3.contains("expansionRatio"))
        lyrOne.expansionRatio = jptch3["expansionRatio"].as<int>();
      else {
        lyrOne.expansionRatio = 1;
      }

      if (jptch3.contains("finalLayerThickness"))
        lyrOne.finalLayerThickness = jptch3["finalLayerThickness"].as<int>();
      else {
        lyrOne.finalLayerThickness = 1;
      }

      if (jptch3.contains("minThickness"))
        lyrOne.minThickness = jptch3["minThickness"].as<int>();
      else {
        lyrOne.minThickness = 1;
      }
      params->layerVec.push_back(lyrOne);
    }
  }

  // Mesh Quality Controls
  if (qcMeshParams.contains("maxNonOrtho"))
    params->qcMaxNOrtho = qcMeshParams["maxNonOrtho"].as<int>();
  else {
    params->qcMaxNOrtho = 65;
  }
  if (qcMeshParams.contains("maxBoundarySkewness"))
    params->qcMaxBndrySkew = qcMeshParams["maxBoundarySkewness"].as<double>();
  else {
    params->qcMaxBndrySkew = 20;
  }
  if (qcMeshParams.contains("maxInternalSkewness"))
    params->qcMaxIntSkew = qcMeshParams["maxInternalSkewness"].as<double>();
  else {
    params->qcMaxIntSkew = 4;
  }
  if (qcMeshParams.contains("maxConcave"))
    params->qcMaxConc = qcMeshParams["maxConcave"].as<double>();
  else {
    params->qcMaxConc = 80;
  }
  if (qcMeshParams.contains("minVol"))
    params->qcMinVol = qcMeshParams["minVol"].as<double>();
  else {
    params->qcMinVol = 1e-13;
  }
  if (qcMeshParams.contains("minTetQuality"))
    params->qcMinTetQ = qcMeshParams["minTetQuality"].as<double>();
  else {
    params->qcMinTetQ = 1e-15;
  }
  if (qcMeshParams.contains("minArea"))
    params->qcMinArea = qcMeshParams["minArea"].as<double>();
  else {
    params->qcMinArea = -1;
  }
  if (qcMeshParams.contains("minTwist"))
    params->qcMinTwist = qcMeshParams["minTwist"].as<double>();
  else {
    params->qcMinTwist = 0.02;
  }
  if (qcMeshParams.contains("minFaceWeight"))
    params->qcMinFaceW = qcMeshParams["minFaceWeight"].as<double>();
  else {
    params->qcMinFaceW = 0.05;
  }
  if (qcMeshParams.contains("minVolRatio"))
    params->qcMinVolRto = qcMeshParams["minVolRatio"].as<double>();
  else {
    params->qcMinVolRto = 0.01;
  }
  if (qcMeshParams.contains("minDeterminant"))
    params->qcMinDet = qcMeshParams["minDeterminant"].as<double>();
  else {
    params->qcMinDet = 0.001;
  }
  if (qcMeshParams.contains("minTriangleTwist"))
    params->qcMinTrTwist = qcMeshParams["minTriangleTwist"].as<double>();
  else {
    params->qcMinTrTwist = -1;
  }
  if (qcMeshParams.contains("qcnSmoothScale"))
    params->qcSmthScale = qcMeshParams["qcnSmoothScale"].as<int>();
  else {
    params->qcSmthScale = 5;
  }
  if (qcMeshParams.contains("errorReduction"))
    params->qcErrRedctn = qcMeshParams["errorReduction"].as<double>();
  else {
    params->qcErrRedctn = 0.75;
  }
  if (shmparams.contains("mergeTolerance"))
    params->mergeTol = shmparams["mergeTolerance"].as<double>();
  else {
    params->mergeTol = 1e-06;
  }

  // Needs attention
  std::string cap2 = "GeomRefinementRegions";
  if (shmparams.contains(cap2)) {
    params->_withGeomRefReg = true;

    for (auto jptch2 : shmparams[cap2].array_range()) {
      shmRegionRefine geomRef;

      if (jptch2.contains("PatchName"))
        geomRef.patchNm = jptch2["PatchName"].as<std::string>();
      else {
        std::cerr << "Please define \"PatchName\"\n" << std::endl;
        throw;
      }
      if (jptch2.contains("Mode"))
        geomRef.mode = jptch2["Mode"].as<std::string>();
      else {
        std::cerr << "Please define \"Mode\"\n" << std::endl;
        throw;
      }
      if (jptch2.contains("MinLevel"))
        geomRef.minLvl = jptch2["MinLevel"].as<int>();
      else {
        std::cerr << "Please define \"MinLevel\"\n" << std::endl;
        throw;
      }
      if (jptch2.contains("MaxLevel"))
        geomRef.maxLvl = jptch2["MaxLevel"].as<int>();
      else {
        std::cerr << "Please define \"MaxLevel\"\n" << std::endl;
        throw;
      }

      (params->geomRefs).push_back(geomRef);
    }
  }

  // blockmesh parameters
  blockMeshParams *bmparams = new blockMeshParams();
  std::string defaults3 =
      inputjson["Meshing Parameters"]["blockMesh Parameters"].as<std::string>();
  jsoncons::json bmshparams =
      inputjson["Meshing Parameters"]["blockMesh Parameters"];

  if (inputjson["Mesh File Options"].contains("Input Dict File"))
    bmparams->_ownBlockMshDict =
        inputjson["Mesh File Options"]["Input Dict File"].as<bool>();
  // parameter parsing starts here
  if (bmshparams.contains("Block Geometry"))
    bmparams->_isBlock = bmshparams["Block Geometry"].as<bool>();
  if (bmshparams.contains("Sphere Geometry"))
    bmparams->_isSphere = bmshparams["Sphere Geometry"].as<bool>();
  if (bmshparams.contains("Cylinder/Tapered_Cone Geometry"))
    bmparams->_isCylinder_TCone =
        bmshparams["Cylinder/Tapered_Cone Geometry"].as<bool>();
  if (bmshparams.contains("scaleToMeters"))
    bmparams->cnvrtToMeters = bmshparams["scaleToMeters"].as<double>();
  if (bmshparams.contains("XdirectionCells"))
    bmparams->cellsXDir = bmshparams["XdirectionCells"].as<int>();
  if (bmshparams.contains("YdirectionCells"))
    bmparams->cellsYDir = bmshparams["YdirectionCells"].as<int>();
  if (bmshparams.contains("ZdirectionCells"))
    bmparams->cellsZDir = bmshparams["ZdirectionCells"].as<int>();

  if (bmshparams.contains("Block Parameters")) {
    bmparams->_autoGenerateBox = false;
    if (bmshparams["Block Parameters"].contains("X1"))
      bmparams->initX = bmshparams["Block Parameters"]["X1"].as<double>();
    if (bmshparams["Block Parameters"].contains("Y1"))
      bmparams->initY = bmshparams["Block Parameters"]["Y1"].as<double>();
    if (bmshparams["Block Parameters"].contains("Z1"))
      bmparams->initZ = bmshparams["Block Parameters"]["Z1"].as<double>();
    if (bmshparams["Block Parameters"].contains("LengthX"))
      bmparams->lenX = bmshparams["Block Parameters"]["LengthX"].as<double>();
    if (bmshparams["Block Parameters"].contains("LengthY"))
      bmparams->lenY = bmshparams["Block Parameters"]["LengthY"].as<double>();
    if (bmshparams["Block Parameters"].contains("LengthZ"))
      bmparams->lenZ = bmshparams["Block Parameters"]["LengthZ"].as<double>();
    if (bmshparams["Block Parameters"].contains("GradingXdir"))
      bmparams->smplGradingX =
          bmshparams["Block Parameters"]["GradingXdir"].as<double>();
    if (bmshparams["Block Parameters"].contains("GradingYdir"))
      bmparams->smplGradingY =
          bmshparams["Block Parameters"]["GradingYdir"].as<double>();
    if (bmshparams["Block Parameters"].contains("GradingZdir"))
      bmparams->smplGradingZ =
          bmshparams["Block Parameters"]["GradingZdir"].as<double>();
  }

  if (bmshparams.contains("Sphere Parameters")) {
    if (bmshparams["Sphere Parameters"].contains("Center X"))
      bmparams->centerX =
          bmshparams["Sphere Parameters"]["Center X"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("Center Y"))
      bmparams->centerY =
          bmshparams["Sphere Parameters"]["Center Y"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("Center Z"))
      bmparams->centerZ =
          bmshparams["Sphere Parameters"]["Center Z"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("Radius"))
      bmparams->radius = bmshparams["Sphere Parameters"]["Radius"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
      bmparams->sphrGradingX =
          bmshparams["Sphere Parameters"]["GradingXdir"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("GradingYdir"))
      bmparams->sphrGradingY =
          bmshparams["Sphere Parameters"]["GradingYdir"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
      bmparams->sphrGradingZ =
          bmshparams["Sphere Parameters"]["GradingZdir"].as<double>();
  }

  if (bmshparams.contains("Cylinder/Tapered_Cone Parameters")) {
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center X"))
      bmparams->centerCyl[0] =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Center X"]
              .as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Y"))
      bmparams->centerCyl[1] =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Y"]
              .as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Z"))
      bmparams->centerCyl[2] =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Z"]
              .as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Radius1"))
      bmparams->radius1 =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"]
              .as<double>();

    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Radius2")) {
      bmparams->radius2 =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius2"]
              .as<double>();
    } else {
      bmparams->radius2 =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"]
              .as<double>();
    }

    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("GradingXdir"))
      bmparams->cylGrading[0] =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingXdir"]
              .as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("GradingYdir"))
      bmparams->cylGrading[1] =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingYdir"]
              .as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("GradingXdir"))
      bmparams->cylGrading[2] =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingZdir"]
              .as<double>();

    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Height"))
      bmparams->height =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Height"].as<double>();
  }

  MeshManipulationFoamParams *mparams = new MeshManipulationFoamParams();
  std::string defaults4 =
      inputjson["Meshing Parameters"]["MeshManipulation Parameters"]
          .as<std::string>();
  jsoncons::json pmshparams =
      inputjson["Meshing Parameters"]["MeshManipulation Parameters"];

  // parameter parsing starts here
  if (pmshparams.contains("Enable SurfLambdaMuSmooth"))
    mparams->_doSurfaceLMSmth =
        pmshparams["Enable SurfLambdaMuSmooth"].as<bool>();
  if (pmshparams.contains("Enable splitMeshRegions"))
    mparams->_doSplitMshRegs = pmshparams["Enable splitMeshRegions"].as<bool>();
  if (pmshparams.contains("Enable MergeMeshes"))
    mparams->_doMergeMsh = pmshparams["Enable MergeMeshes"].as<bool>();
  if (pmshparams.contains("Enable CreatePatch"))
    mparams->_doCreatePtchs = pmshparams["Enable CreatePatch"].as<bool>();
  if (pmshparams.contains("Enable foamToSurface"))
    mparams->_doFoam2Surf = pmshparams["Enable foamToSurface"].as<bool>();
  if (pmshparams.contains("Enable surfaceSplitByTopology"))
    mparams->_doSurfSplit =
        pmshparams["Enable surfaceSplitByTopology"].as<bool>();

  if (pmshparams.contains("SurfLambdaMuSmooth Parameters")) {
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("AddFeatureFile?"))
      mparams->_addFeatureFile =
          pmshparams["SurfLambdaMuSmooth Parameters"]["AddFeatureFile?"]
              .as<bool>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Input STL File"))
      mparams->slmssurfaceFile =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Input STL File"]
              .as<std::string>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Output STL File"))
      mparams->slmsoutputFile =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Output STL File"]
              .as<std::string>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Lambda"))
      mparams->lambda =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Lambda"].as<double>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Mu"))
      mparams->mu =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Mu"].as<double>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains(
            "Smoothing Interations"))
      mparams->slmsIterations =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Smoothing Interations"]
              .as<int>();
  }

  if (pmshparams.contains("splitMeshRegions Parameters")) {
    if (pmshparams["splitMeshRegions Parameters"].contains("overwrite?"))
      mparams->_overwriteMsh =
          pmshparams["splitMeshRegions Parameters"]["overwrite?"].as<bool>();
    if (pmshparams["splitMeshRegions Parameters"].contains("usecellZones?"))
      mparams->_cellZones =
          pmshparams["splitMeshRegions Parameters"]["usecellZones?"].as<bool>();
  }

  if (pmshparams.contains("mergeMeshes Parameters")) {
    if (pmshparams["mergeMeshes Parameters"].contains("Master Region"))
      mparams->masterCase =
          pmshparams["mergeMeshes Parameters"]["Master Region"]
              .as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("Add Region"))
      mparams->addCase =
          pmshparams["mergeMeshes Parameters"]["Add Region"].as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("overwrite?"))
      mparams->_overwriteMergeMsh =
          pmshparams["mergeMeshes Parameters"]["overwrite?"].as<bool>();
    if (pmshparams["mergeMeshes Parameters"].contains("Master Region Path"))
      mparams->masterCasePath =
          pmshparams["mergeMeshes Parameters"]["Master Region Path"]
              .as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("Add Region Path"))
      mparams->addCasePath =
          pmshparams["mergeMeshes Parameters"]["Add Region Path"]
              .as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("Number of Domains"))
      mparams->numDomains =
          pmshparams["mergeMeshes Parameters"]["Number of Domains"].as<int>();
  }

  if (pmshparams.contains("createPatch Parameters")) {
    if (pmshparams["createPatch Parameters"].contains("Surrounding PatchName"))
      mparams->surroundingName =
          pmshparams["createPatch Parameters"]["Surrounding PatchName"]
              .as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Packs PatchName"))
      mparams->packsName =
          pmshparams["createPatch Parameters"]["Packs PatchName"]
              .as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Surrounding PatchType"))
      mparams->srrndngPatchType =
          pmshparams["createPatch Parameters"]["Surrounding PatchType"]
              .as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Packs PatchType"))
      mparams->packsPatchType =
          pmshparams["createPatch Parameters"]["Packs PatchType"]
              .as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("overwrite?"))
      mparams->_overwritecpMsh =
          pmshparams["createPatch Parameters"]["overwrite?"].as<bool>();
    if (pmshparams["createPatch Parameters"].contains("Packs Path"))
      mparams->pathPacks =
          pmshparams["createPatch Parameters"]["Packs Path"].as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Surrounding Path"))
      mparams->pathSurrounding =
          pmshparams["createPatch Parameters"]["Surrounding Path"]
              .as<std::string>();
  }

  if (pmshparams.contains("foamToSurface Parameters")) {
    if (pmshparams["foamToSurface Parameters"].contains("Output File Path"))
      mparams->outSurfName =
          pmshparams["foamToSurface Parameters"]["Output File Path"]
              .as<std::string>();
  }

  if (pmshparams.contains("surfaceSplitByTopology Parameters")) {
    if (pmshparams["surfaceSplitByTopology Parameters"].contains("Input File"))
      mparams->surfFile =
          pmshparams["surfaceSplitByTopology Parameters"]["Input File"]
              .as<std::string>();
    if (pmshparams["surfaceSplitByTopology Parameters"].contains("Output File"))
      mparams->outSurfFile =
          pmshparams["surfaceSplitByTopology Parameters"]["Output File"]
              .as<std::string>();
  }

  // creating triSurface directory for workflow
  std::string dir_path = "./constant";
  boost::filesystem::path dir(dir_path);
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr
        << "Problem in creating triSurface directory for the snappyHexMesh"
        << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  // directory initialization
  dir = "./constant/triSurface";
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr
        << "Problem in creating triSurface directory for the snappyHexMesh"
        << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  MeshManipulationFoam *objMsh = new MeshManipulationFoam(mparams);
  const char *nameFile = "a";

  /*
  // only modified version of the cfmesh can be rolled in
  // cfmesh
  cfmeshGen* objCF = new cfmeshGen(cfparams);
  objCF->createMeshFromSTL(nameFile);

  // foamToSurface
  objMsh->foamToSurface();
  */

  int nDom = 2;

  // blockMesh
  blockMeshGen *objBM = new blockMeshGen(bmparams);
  objBM->createMeshFromSTL(nameFile);

  // snappyHexMesh
  snappymeshGen *objSHM = new snappymeshGen(params);
  objSHM->createMeshFromSTL(nameFile);

  // splitMeshRegions
  std::pair<int, int> dirStat =
      objMsh->splitMshRegions();  // outputs the region number
                                  // skipped during splitting process
  int skippedDir = dirStat.first;
  int totalRegs = dirStat.second - 1;

  // mergeMeshes
  objMsh->mergeMeshes(skippedDir, nDom);

  // createPatch
  objMsh->createPatch(skippedDir);

  // read current mesh and write it to separate VTK/VTU files
  bool readDB = false;
  // converts pack mesh
  std::string regNme;
  if (skippedDir == 1)
    regNme = "domain2";
  else
    regNme = "domain1";

  meshBase *fm = new FOAM::foamMesh(readDB);
  fm->read(regNme);
  vtkMesh *vm = new vtkMesh(fm->getDataSet(), ofname1);
  vm->report();
  vm->write();

  // converts surronding mesh
  regNme = "domain0";
  meshBase *fm2 = new FOAM::foamMesh(readDB);
  fm2->read(regNme);
  vtkMesh *vm2 = new vtkMesh(fm2->getDataSet(), ofname2);
  vm2->report();
  vm2->write();

  // cleaning up
  if (vm) delete vm;
  if (fm) delete fm;
  if (vm2) delete vm2;
  if (fm2) delete fm2;
  if (objMsh) delete objMsh;
  if (objSHM) delete objSHM;
  // if (objCF)
  //  delete objCF;
  if (objBM) delete objBM;

  return 0;
}

// TEST macros
TEST(PackMeshing, Generation) { EXPECT_EQ(0, generate(inp_json)); }

TEST(PackMeshing, NumberOfNodesnCellsPacks) {
  if (ref) delete ref;
  ref = meshBase::Create(inputjson["Pack Reference File"].as<std::string>());
  meshBase *cmp1 = meshBase::Create("geom_pack_mesh.vtu");
  EXPECT_TRUE((cmp1->getNumberOfPoints() == ref->getNumberOfPoints()) &&
              (cmp1->getNumberOfCells() == ref->getNumberOfCells()));
}

TEST(PackMeshing, NumberOfNodesNodesnCellsSurrounding) {
  if (ref) delete ref;
  ref = meshBase::Create(
      inputjson["Surrounding Reference File"].as<std::string>());
  meshBase *cmp2 = meshBase::Create("geom_surrounding_mesh.vtu");
  EXPECT_TRUE((cmp2->getNumberOfPoints() == ref->getNumberOfPoints()) &&
              (cmp2->getNumberOfCells() == ref->getNumberOfCells()));
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 1);
  inp_json = argv[1];

  if (!inp_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }

  // running tests
  int res = RUN_ALL_TESTS();

  // clean up
  if (mesh) delete mesh;

  return res;
}
