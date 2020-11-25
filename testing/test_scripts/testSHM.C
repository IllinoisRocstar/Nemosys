#include <gtest.h>
#include <algorithm>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iterator>
#include <string>
#include "NemDriver.H"
#include "snappymeshGen.H"
#include "snappymeshParams.H"
#include "vtkMesh.H"
namespace fs = boost::filesystem;

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

void copyDirectoryRecursively(const fs::path &sourceDir,
                              const fs::path &destinationDir) {
  if (!fs::exists(sourceDir) || !fs::is_directory(sourceDir)) {
    throw std::runtime_error("Source directory " + sourceDir.string() +
                             " does not exist or is not a directory");
  }
  if (fs::exists(destinationDir)) {
    throw std::runtime_error("Destination directory " +
                             destinationDir.string() + " already exists");
  }
  if (!fs::create_directory(destinationDir)) {
    throw std::runtime_error("Cannot create destination directory " +
                             destinationDir.string());
  }

  for (const auto &dirEnt : fs::recursive_directory_iterator{sourceDir}) {
    const auto &path = dirEnt.path();
    auto relativePathStr = path.string();
    boost::algorithm::replace_first(relativePathStr, sourceDir.string(), "");
    fs::copy(path, destinationDir / relativePathStr);
  }
}

void snappyLogistics() {
  const char dir_path[] = "./constant/polyMesh";
  const char cpydir_path[] = "./constant/polyMesh.Orig";
  boost::filesystem::path dir1(dir_path);
  boost::filesystem::path dir2(cpydir_path);
  boost::filesystem::remove_all(dir1);
  copyDirectoryRecursively(dir2, dir1);

  return;
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

  snappymeshParams *params = new snappymeshParams();
  std::string ifname =
      inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  std::string ofname =
      inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  jsoncons::json shmparams =
      inputjson["Meshing Parameters"]["snappyHexMesh Parameters"];

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
    // params->_withSurfRefReg = true;

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
        surfRef.minLvl = jptch3["MinLevel"].as<double>();
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

  // Refinement Regions
  std::string cap2 = "GeomRefinementRegions";
  if (shmparams.contains(cap2)) {
    // params->_withGeomRefReg = true;

    for (auto jptch2 : shmparams[cap2].array_range()) {
      shmRegionRefine geomRef;

      if (jptch2.contains("Patch Name"))
        geomRef.patchNm = jptch2["Patch Name"].as<std::string>();
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
        geomRef.minLvl = jptch2["MinLevel"].as<double>();
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
        ftrOne.minLvl = jptch3["MinLevel"].as<double>();
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
  if (layerParams.contains("firstLayerThickness"))
    params->firstLyrThickness = layerParams["firstLayerThickness"].as<double>();
  else
    params->firstLyrThickness = -1.0;
  if (layerParams.contains("thickness"))
    params->thickness = layerParams["thickness"].as<double>();
  else
    params->thickness = -1.0;
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
  if (layerParams.contains("nMedialAxisIter"))
    params->nMedialAxisIter = layerParams["nMedialAxisIter"].as<int>();
  else
    params->nMedialAxisIter = -1;
  if (layerParams.contains("nSmoothDisplacement"))
    params->nSmoothDisplacement = layerParams["nSmoothDisplacement"].as<int>();
  else
    params->nSmoothDisplacement = -1;

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

      if (jptch3.contains("firstLayerThickness"))
        lyrOne.firstLyrThickness = jptch3["firstLayerThickness"].as<double>();
      else {
        lyrOne.firstLyrThickness = -1.0;
      }

      if (jptch3.contains("thickness"))
        lyrOne.finalLayerThickness = jptch3["thickness"].as<double>();
      else {
        lyrOne.thickness = -1.0;
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

  snappymeshGen *generator =
      new snappymeshGen(dynamic_cast<snappymeshParams *>(params));
  generator->createMeshFromSTL("");
  mesh = vtkMesh::Create(generator->getDataSet(), ofname);
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();

  return 0;
}

// Write a new name for snappyHexMesh
// TEST macros
TEST(snappyHexMesh, Generation) { EXPECT_EQ(0, generate(inp_json)); }

TEST(snappyHexMesh, NumberOfNodes) {
  if (ref) delete ref;
  ref = meshBase::Create(inputjson["Reference File"].as<std::string>());
  EXPECT_EQ(mesh->getNumberOfPoints(), ref->getNumberOfPoints());
}

TEST(snappyHexMesh, NumberOfCells) {
  if (ref) delete ref;
  ref = meshBase::Create(inputjson["Reference File"].as<std::string>());
  EXPECT_EQ(mesh->getNumberOfCells(), ref->getNumberOfCells());
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 1);
  inp_json = argv[1];

  if (!inp_json) {
    std::cerr << "No input file defined" << std::endl;
  }

  int res = RUN_ALL_TESTS();

  // clean up
  if (mesh) delete mesh;

  snappyLogistics();  // clean up and replace mesh

  return res;
}
