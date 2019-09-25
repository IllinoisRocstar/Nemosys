#include "NemDriver.H"
#include "cfmeshGen.H"
#include "cfmeshParams.H"
#include "blockMeshGen.H"
#include "blockMeshParams.H"
#include "MeshManipulationFoam.H"
#include "MeshManipulationFoamParams.H"
#include "snappymeshGen.H"
#include "snappymeshParams.H"
#include <PackMeshDriver.H>
#include "vtkMesh.H"
#include <gtest.h>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/filesystem.hpp>

const char* inp_json;
meshBase* mesh;
meshBase* ref;
jsoncons::json inputjson;

// Aux functions
bool compareFiles(const std::string& p1, const std::string& p2) {
  std::ifstream f1(p1, std::ifstream::binary|std::ifstream::ate);
  std::ifstream f2(p2, std::ifstream::binary|std::ifstream::ate);

  if (f1.fail() || f2.fail()) {
    return false; //file problem
  }

  if (f1.tellg() != f2.tellg()) {
    return false; //size mismatch
  }

  //seek back to beginning and use std::equal to compare contents
  f1.seekg(0, std::ifstream::beg);
  f2.seekg(0, std::ifstream::beg);
  return std::equal(std::istreambuf_iterator<char>(f1.rdbuf()),
                    std::istreambuf_iterator<char>(),
                    std::istreambuf_iterator<char>(f2.rdbuf()));
}

// Test implementations
int generate(const char* jsonF)
{
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if(!inputStream.good())
  {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }
  
  jsoncons::json inputjson_tmp;
  inputStream >> inputjson_tmp;
  inputjson = inputjson_tmp[0]; 

  // TODO: refactor the input processing, use actual driver objtect
  //       intead of parsing here.
  std::string ifname = inputjson["Mesh File Options"]
                                ["Input Geometry File"].as<std::string>();
  std::string ofname1 = inputjson["Mesh File Options"]
                                ["Output Pack Mesh File"].as<std::string>();
  std::string ofname2 = inputjson["Mesh File Options"]
                                ["Output Surrounding Mesh File"].as<std::string>();

  cfmeshParams* cfparams = new cfmeshParams();
  std::string defaults1 = 
      inputjson["Meshing Parameters"]["CFMesh Parameters"].as<std::string>();
  jsoncons::json cfmparams = inputjson["Meshing Parameters"]["CFMesh Parameters"];

  // required params here
  // cad file
  if (inputjson["Mesh File Options"].contains("Input Geometry File"))
      cfparams->geomFilePath = 
          inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else
  {
      std::cerr << "A geometry file should be supplied.\n";
      throw;        
  }
  
  // mesh generator
  if (cfmparams.contains("Generator"))
      cfparams->generator = cfmparams["Generator"].as<std::string>();
  else
  {
      std::cerr << "A mesh generation method should be selected.\n";
      std::cerr << "Options: cartesian2D tetMesh\n";
      throw;        
  }
  
  // rest of params are optional
  if (cfmparams.contains("MaxCellSize"))
    cfparams->maxCellSize = 
        cfmparams["MaxCellSize"].as<double>();
  if (cfmparams.contains("MinCellSize"))
    cfparams->minCellSize = 
        cfmparams["MinCellSize"].as<double>();
  if (cfmparams.contains("BoundaryCellSize"))
    cfparams->bndryCellSize = 
        cfmparams["BoundaryCellSize"].as<double>();
  if (cfmparams.contains("KeepCellsIntersectingBoundary"))
    cfparams->keepCellIB = 
        cfmparams["KeepCellsIntersectingBoundary"].as<double>();
  if (cfmparams.contains("CheckForGluedMesh"))
    cfparams->chkGluMsh = 
        cfmparams["CheckForGluedMesh"].as<double>();
  if (cfmparams.contains("AllowDisconnectedDomains"))
      cfparams->_alwDiscDomains = 
        cfmparams["AllowDisconnectedDomains"].as<bool>();


  // optional capability
  std::string cap = "BoundaryLayers";
  if (cfmparams.contains(cap))
  {
      cfparams->_withBndLyr = true;
      cfparams->blNLyr = cfmparams[cap]["NLayers"].as<double>();
      cfparams->blThkRto = cfmparams[cap]["ThicknessRatio"].as<double>();
      if (cfmparams[cap].contains("MaxFirstLayerThickness"))
          cfparams->maxFrstLyrThk = cfmparams[cap]["MaxFirstLayerThickness"].as<double>();
      if (cfmparams[cap].contains("AllowDiscontinuity"))
          cfparams->alwDiscont = cfmparams[cap]["AllowDiscontinuity"].as<bool>();

      // patch boundary layers
      std::string subcap = "PatchBoundaryLayers";
      if (cfmparams[cap].contains(subcap))
      {
        cfparams->_withBndLyrPtch = true;
        for (auto jptch : cfmparams[cap][subcap].array_range())
        {
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
  if (cfmparams.contains(cap))
  {
      cfparams->_withSrfEdg = true;
      cfparams->srfEdgAng = cfmparams[cap]["Angle"].as<double>();
  }

  // optional capability
  cap = "ObjectRefinements";
  if (cfmparams.contains(cap))
  {
    cfparams->_withObjRfn = true;
    for (auto refObj : cfmparams[cap].array_range())
    {
      cfmObjRef objRef;
      objRef.name= refObj["Name"].as<std::string>();
      for (const auto&  prm: refObj["Params"].object_range())
      {
          std::string key = std::string(prm.key());
          std::string val = prm.value().as<std::string>();
          objRef.params[key] = val;
      }
      (cfparams->objRefLst).push_back(objRef);
    }
  }

  // optional capability
  cap = "ImproveMeshQuality";
  if (cfmparams.contains(cap))
  {
    cfparams->_withMshQlt = true;
    cfparams->qltNItr = cfmparams[cap]["NIterations"].as<int>();
    cfparams->qltNLop = cfmparams[cap]["NLoops"].as<int>();
    cfparams->qltQltThr = cfmparams[cap]["QualityThreshold"].as<double>();
    cfparams->qltNSrfItr = cfmparams[cap]["NSurfaceIterations"].as<int>();
    cfparams->qltConCelSet = 
      cfmparams[cap].get_with_default("ConstrainedCellsSet","none");
  }

  // optional capability
  cap = "LocalRefinement";
  if (cfmparams.contains(cap))
  {
    cfparams->_withLclRef = true;
    for (auto jptch : cfmparams[cap].array_range())
    {
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
  if (cfmparams.contains(cap))
  {
    cfparams->_withRenBndry = true;
    cfmRenBndry renBndry;
    renBndry.defName = cfmparams[cap]["DefaultName"].as<std::string>();
    renBndry.defType = cfmparams[cap]["DefaultType"].as<std::string>();
    for (auto jnw : cfmparams[cap]["NewPatchNames"].array_range())
    {
      cfmNewPatch nwPatch = std::make_tuple( 
          jnw["Name"].as<std::string>(),
          jnw["NewName"].as<std::string>(),
          jnw["NewType"].as<std::string>() );
      renBndry.newPatches.push_back(nwPatch);
    }
    (cfparams->renBndry) = renBndry;
  }


  // snappy
  snappymeshParams* snappyparams = new snappymeshParams();
  std::string defaults2 = 
  inputjson["Meshing Parameters"]["snappyHexMesh Parameters"].as<std::string>();
  jsoncons::json shmparams = inputjson["Meshing Parameters"]["snappyHexMesh Parameters"];
  if (inputjson["Mesh File Options"].contains("Input Geometry File"))
      snappyparams->geomFileName = 
          inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else
  {
    std::cerr << "A geometry file should be supplied.\n";
    throw;        
  }

  if (shmparams.contains("InputPatchName"))
      snappyparams->geomPatchName = shmparams["InputPatchName"].as<std::string>();
  else
  {
    std::cerr << "A patch name for input geometry must be defined.\n";
    throw;        
  }

  if (shmparams.contains("SurfPatchName"))
      snappyparams->surfRefPatch = shmparams["SurfPatchName"].as<std::string>();
  else
  {
    std::cerr << "A patch name for surface refinement must be defined.\n";
    throw;        
  }
  
  if (shmparams.contains("Castellated Mesh"))
    snappyparams->_withCastMesh =
      shmparams["Castellated Mesh"].as<bool>();
  if (shmparams.contains("Snapping"))
    snappyparams->_withSnap =
      shmparams["Snapping"].as<bool>();
  if (shmparams.contains("Layer Addition"))
    snappyparams->_withLayers =
      shmparams["Layer Addition"].as<bool>();
  if (shmparams.contains("CellZones"))
    snappyparams->_withCellZones =
      shmparams["CellZones"].as<bool>();
  if (shmparams.contains("RegionRefine"))
    snappyparams->_withGeomRefReg =
      shmparams["RegionRefine"].as<bool>();
  if (shmparams.contains("SurfaceRefine"))
    snappyparams->_withSurfRefReg =
      shmparams["SurfaceRefine"].as<bool>();
  if (shmparams.contains("maxLocalCells"))
    snappyparams->maxLCells =
      shmparams["maxLocalCells"].as<int>();
  if (shmparams.contains("maxGlobalCells"))
    snappyparams->maxGCells =
      shmparams["maxGlobalCells"].as<int>();
  if (shmparams.contains("minRefCells"))
    snappyparams->minRefCells =
      shmparams["minRefCells"].as<int>();
  if (shmparams.contains("nCellsBetweenLevels"))
    snappyparams->cellsBetnLvls =
      shmparams["nCellsBetweenLevels"].as<int>();
  if (shmparams.contains("surfaceRefinementLvlMin"))
    snappyparams->refSurfLvlMin =
      shmparams["surfaceRefinementLvlMin"].as<int>();
  if (shmparams.contains("surfaceRefinementLvlMax"))
    snappyparams->refSurfLvlMax =
      shmparams["surfaceRefinementLvlMax"].as<int>();
  if (shmparams.contains("resolveFeatureAngle"))
    snappyparams->featAngle =
      shmparams["resolveFeatureAngle"].as<double>();
  if (shmparams.contains("locationInMeshX"))
    snappyparams->locMeshX =
      shmparams["locationInMeshX"].as<double>();
  if (shmparams.contains("locationInMeshY"))
    snappyparams->locMeshY =
      shmparams["locationInMeshY"].as<double>();
  if (shmparams.contains("locationInMeshZ"))
    snappyparams->locMeshZ =
      shmparams["locationInMeshZ"].as<double>();
  if (shmparams.contains("allowFreeStandingZoneFaces"))
    snappyparams->_alwFreeZone =
      shmparams["allowFreeStandingZoneFaces"].as<double>();
  if (shmparams.contains("nSmoothPatch"))
    snappyparams->snapSmthPatch =
      shmparams["nSmoothPatch"].as<int>();
  if (shmparams.contains("tolerance"))
    snappyparams->snapTol =
      shmparams["tolerance"].as<double>();
  if (shmparams.contains("snapSolveIter"))
    snappyparams->solveSnapIter =
      shmparams["snapSolveIter"].as<int>();
  if (shmparams.contains("snapRelaxIter"))
    snappyparams->relaxSnapIter =
      shmparams["snapRelaxIter"].as<int>();
  if (shmparams.contains("relativeSizes"))
    snappyparams->_relSize =
      shmparams["relativeSizes"].as<double>();
  if (shmparams.contains("expansionRatio"))
    snappyparams->expRatio =
      shmparams["expansionRatio"].as<double>();
  if (shmparams.contains("finalLayerThickness"))
    snappyparams->finLThick =
      shmparams["finalLayerThickness"].as<double>();
  if (shmparams.contains("minThickness"))
    snappyparams->minThick =
      shmparams["minThickness"].as<double>();
  if (shmparams.contains("nGrow"))
    snappyparams->nGrow =
      shmparams["nGrow"].as<int>();
  if (shmparams.contains("featureAngle"))
    snappyparams->lyrFeatAngle =
      shmparams["featureAngle"].as<double>();
  if (shmparams.contains("nRelaxIter"))
    snappyparams->lyrRelaxIter =
      shmparams["nRelaxIter"].as<int>();
  if (shmparams.contains("nSmoothSurfaceNormals"))
    snappyparams->lyrSmthSurfNorm =
      shmparams["nSmoothSurfaceNormals"].as<int>();
  if (shmparams.contains("nSmoothNormals"))
    snappyparams->lyrSmthNorm =
      shmparams["nSmoothNormals"].as<int>();
  if (shmparams.contains("nSmoothThickness"))
    snappyparams->lyrSmthThick =
      shmparams["nSmoothThickness"].as<int>();
  if (shmparams.contains("maxFaceThicknessRatio"))
    snappyparams->lyrMaxFcTR =
      shmparams["maxFaceThicknessRatio"].as<double>();
  if (shmparams.contains("maxThicknessToMedialRatio"))
    snappyparams->lyrMaxThickTMR =
      shmparams["maxThicknessToMedialRatio"].as<double>();
  if (shmparams.contains("minMedialAxisAngle"))
    snappyparams->lyrMinMedAngl =
      shmparams["minMedialAxisAngle"].as<double>();
  if (shmparams.contains("nBufferCellsNoExtrude"))
    snappyparams->lyrBuffrCells =
      shmparams["nBufferCellsNoExtrude"].as<int>();
  if (shmparams.contains("nLayerIter"))
    snappyparams->lyrIter =
      shmparams["nLayerIter"].as<int>();
  if (shmparams.contains("maxNonOrtho"))
    snappyparams->qcMaxNOrtho =
      shmparams["maxNonOrtho"].as<int>();
  if (shmparams.contains("maxBoundarySkewness"))
    snappyparams->qcMaxBndrySkew =
      shmparams["maxBoundarySkewness"].as<double>();
  if (shmparams.contains("maxInternalSkewness"))
    snappyparams->qcMaxIntSkew =
      shmparams["maxInternalSkewness"].as<double>();
  if (shmparams.contains("maxConcave"))
    snappyparams->qcMaxConc =
      shmparams["maxConcave"].as<double>();
  if (shmparams.contains("minVol"))
    snappyparams->qcMinVol =
      shmparams["minVol"].as<double>();
  if (shmparams.contains("minTetQuality"))
    snappyparams->qcMinTetQ =
      shmparams["minTetQuality"].as<double>();
  if (shmparams.contains("minArea"))
    snappyparams->qcMinArea =
      shmparams["minArea"].as<double>();
  if (shmparams.contains("minTwist"))
    snappyparams->qcMinTwist =
      shmparams["minTwist"].as<double>();
  if (shmparams.contains("minFaceWeight"))
    snappyparams->qcMinFaceW =
      shmparams["minFaceWeight"].as<double>();
  if (shmparams.contains("minVolRatio"))
    snappyparams->qcMinVolRto =
      shmparams["minVolRatio"].as<double>();
  if (shmparams.contains("minDeterminant"))
    snappyparams->qcMinDet =
      shmparams["minDeterminant"].as<double>();
  if (shmparams.contains("minTriangleTwist"))
    snappyparams->qcMinTrTwist =
      shmparams["minTriangleTwist"].as<double>();
  if (shmparams.contains("qcnSmoothScale"))
    snappyparams->qcSmthScale =
      shmparams["qcnSmoothScale"].as<int>();
  if (shmparams.contains("errorReduction"))
    snappyparams->qcErrRedctn =
      shmparams["errorReduction"].as<double>();
  if (shmparams.contains("mergeTolerance"))
    snappyparams->mergeTol =
      shmparams["mergeTolerance"].as<double>();

  // optional capability
  std::string cap2 = "GeomRefinementRegions";
  if (shmparams.contains(cap2))
  {
    snappyparams->_withGeomRefReg = true;

    for (auto jptch2 : shmparams[cap2].array_range())
    {
      shmGeomRefine geomRef;

      if (jptch2.contains("PatchName"))
        geomRef.patchNm = jptch2["PatchName"].as<std::string>();

      if (jptch2.contains("searchableShape"))
        geomRef.searchableName = jptch2["searchableShape"].as<std::string>();

      if (jptch2.contains("shapeParams1"))
        geomRef.shapeParameters1 = jptch2["shapeParams1"].as<std::string>();

      if (jptch2.contains("shapeParams2"))
        geomRef.shapeParameters2 = jptch2["shapeParams2"].as<std::string>();

      if (jptch2.contains("Radius"))
        geomRef.rad = jptch2["Radius"].as<double>();

      if (jptch2.contains("Mode"))
        geomRef.mode = jptch2["Mode"].as<std::string>();

      if (jptch2.contains("MinLevel"))
        geomRef.minLvl = jptch2["MinLevel"].as<int>();

      if (jptch2.contains("MaxLevel"))
        geomRef.maxLvl = jptch2["MaxLevel"].as<int>();

      (snappyparams->geomRefs).push_back(geomRef);

    }
  }

  // optional capability
  std::string cap3 = "SurfaceRefinementRegions";
  if (shmparams.contains(cap3))
  {
    snappyparams->_withSurfRefReg = true;
    for (auto jptch3 : shmparams[cap3].array_range())
    {
      shmRegionRef surfRef;
      if (jptch3.contains("PatchName"))
        surfRef.refPatchNm = jptch3["PatchName"].as<std::string>();
      if (jptch3.contains("MinLevel"))
        surfRef.minLvl = jptch3["MinLevel"].as<int>();
      if (jptch3.contains("MaxLevel"))
        surfRef.maxLvl = jptch3["MaxLevel"].as<int>();
      (snappyparams->surfRefs).push_back(surfRef);
    }
  }


  // blockmesh parameters
  blockMeshParams* bmparams = new blockMeshParams();
  std::string defaults3 = 
      inputjson["Meshing Parameters"]["blockMesh Parameters"].as<std::string>();
    jsoncons::json bmshparams = inputjson["Meshing Parameters"]["blockMesh Parameters"];
  
  if (inputjson["Mesh File Options"].contains("Input Dict File"))
    bmparams->_ownBlockMshDict = 
      inputjson["Mesh File Options"]["Input Dict File"].as<bool>();
  // parameter parsing starts here
  if (bmshparams.contains("Block Geometry"))
    bmparams->_isBlock =
      bmshparams["Block Geometry"].as<bool>();
  if (bmshparams.contains("Sphere Geometry"))
    bmparams->_isSphere =
      bmshparams["Sphere Geometry"].as<bool>();
  if (bmshparams.contains("Cylinder/Tapered_Cone Geometry"))
    bmparams->_isCylinder_TCone =
      bmshparams["Cylinder/Tapered_Cone Geometry"].as<bool>();
  if (bmshparams.contains("scaleToMeters"))
    bmparams->cnvrtToMeters = 
      bmshparams["scaleToMeters"].as<double>();
  if (bmshparams.contains("XdirectionCells"))
    bmparams->cellsXDir = 
      bmshparams["XdirectionCells"].as<int>();
  if (bmshparams.contains("YdirectionCells"))
    bmparams->cellsYDir = 
      bmshparams["YdirectionCells"].as<int>();
  if (bmshparams.contains("ZdirectionCells"))
    bmparams->cellsZDir = 
      bmshparams["ZdirectionCells"].as<int>();
    
  if (bmshparams.contains("Block Parameters"))
  {
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
  
  if (bmshparams.contains("Sphere Parameters"))
  {
    if (bmshparams["Sphere Parameters"].contains("Center X"))
      bmparams->centerX = bmshparams["Sphere Parameters"]["Center X"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("Center Y"))
      bmparams->centerY = bmshparams["Sphere Parameters"]["Center Y"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("Center Z"))
      bmparams->centerZ = bmshparams["Sphere Parameters"]["Center Z"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("Radius"))
      bmparams->radius = bmshparams["Sphere Parameters"]["Radius"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
      bmparams->sphrGradingX = bmshparams["Sphere Parameters"]["GradingXdir"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("GradingYdir"))
      bmparams->sphrGradingY = bmshparams["Sphere Parameters"]["GradingYdir"].as<double>();
    if (bmshparams["Sphere Parameters"].contains("GradingXdir"))
      bmparams->sphrGradingZ = bmshparams["Sphere Parameters"]["GradingZdir"].as<double>();
  }
  
  if (bmshparams.contains("Cylinder/Tapered_Cone Parameters"))
  {
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center X"))
      bmparams->centerCyl[0] = bmshparams["Cylinder/Tapered_Cone Parameters"]["Center X"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Y"))
      bmparams->centerCyl[1] = bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Y"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Center Z"))
      bmparams->centerCyl[2] = bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Z"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Radius1"))
      bmparams->radius1 = bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"].as<double>();
    
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Radius2")){
      bmparams->radius2 = bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius2"].as<double>();
    }
    else{
      bmparams->radius2 = bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"].as<double>();
    }
    
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("GradingXdir"))
      bmparams->cylGrading[0] = bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingXdir"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("GradingYdir"))
      bmparams->cylGrading[1] = bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingYdir"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("GradingXdir"))
      bmparams->cylGrading[2] = bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingZdir"].as<double>();
    
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].contains("Height"))
      bmparams->height = bmshparams["Cylinder/Tapered_Cone Parameters"]["Height"].as<double>();
  }


  MeshManipulationFoamParams* mparams = new MeshManipulationFoamParams();
  std::string defaults4 = 
      inputjson["Meshing Parameters"]["MeshManipulation Parameters"].as<std::string>();
    jsoncons::json pmshparams = inputjson["Meshing Parameters"]["MeshManipulation Parameters"];

  
  // parameter parsing starts here
  if (pmshparams.contains("Enable SurfLambdaMuSmooth"))
        mparams->_doSurfaceLMSmth =
          pmshparams["Enable SurfLambdaMuSmooth"].as<bool>();
  if (pmshparams.contains("Enable splitMeshRegions"))
        mparams->_doSplitMshRegs =
          pmshparams["Enable splitMeshRegions"].as<bool>();
  if (pmshparams.contains("Enable MergeMeshes"))
        mparams->_doMergeMsh = 
          pmshparams["Enable MergeMeshes"].as<bool>();
  if (pmshparams.contains("Enable CreatePatch"))
        mparams->_doCreatePtchs = 
          pmshparams["Enable CreatePatch"].as<bool>();
  if (pmshparams.contains("Enable foamToSurface"))
        mparams->_doFoam2Surf = 
          pmshparams["Enable foamToSurface"].as<bool>();
  if (pmshparams.contains("Enable surfaceSplitByTopology"))
        mparams->_doSurfSplit = 
          pmshparams["Enable surfaceSplitByTopology"].as<bool>();

  if (pmshparams.contains("SurfLambdaMuSmooth Parameters"))
  {
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("AddFeatureFile?"))
        mparams->_addFeatureFile =
          pmshparams["SurfLambdaMuSmooth Parameters"]["AddFeatureFile?"].as<bool>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Input STL File"))
        mparams->slmssurfaceFile =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Input STL File"].as<std::string>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Output STL File"))
        mparams->slmsoutputFile =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Output STL File"].as<std::string>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Lambda"))
        mparams->lambda =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Lambda"].as<double>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Mu"))
        mparams->mu =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Mu"].as<double>();
    if (pmshparams["SurfLambdaMuSmooth Parameters"].contains("Smoothing Interations"))
        mparams->slmsIterations =
          pmshparams["SurfLambdaMuSmooth Parameters"]["Smoothing Interations"].as<int>();           
  }


  if (pmshparams.contains("splitMeshRegions Parameters"))
  {
    if (pmshparams["splitMeshRegions Parameters"].contains("overwrite?"))
      mparams->_overwriteMsh = 
        pmshparams["splitMeshRegions Parameters"]["overwrite?"].as<bool>();
    if (pmshparams["splitMeshRegions Parameters"].contains("usecellZones?"))
      mparams->_cellZones = 
        pmshparams["splitMeshRegions Parameters"]["usecellZones?"].as<bool>();
  }


  if (pmshparams.contains("mergeMeshes Parameters"))
  {
    if (pmshparams["mergeMeshes Parameters"].contains("Master Region"))
      mparams->masterCase =
        pmshparams["mergeMeshes Parameters"]["Master Region"].as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("Add Region"))
      mparams->addCase = 
        pmshparams["mergeMeshes Parameters"]["Add Region"].as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("overwrite?"))
      mparams->_overwriteMergeMsh =
        pmshparams["mergeMeshes Parameters"]["overwrite?"].as<bool>();
    if (pmshparams["mergeMeshes Parameters"].contains("Master Region Path"))
      mparams->masterCasePath =
        pmshparams["mergeMeshes Parameters"]["Master Region Path"].as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("Add Region Path"))
      mparams->addCasePath = 
        pmshparams["mergeMeshes Parameters"]["Add Region Path"].as<std::string>();
    if (pmshparams["mergeMeshes Parameters"].contains("Number of Domains"))
      mparams->numDomains = 
        pmshparams["mergeMeshes Parameters"]["Number of Domains"].as<int>();
  }


  if (pmshparams.contains("createPatch Parameters"))
  {
    if (pmshparams["createPatch Parameters"].contains("Surrounding PatchName"))
      mparams->surroundingName = 
        pmshparams["createPatch Parameters"]["Surrounding PatchName"].as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Packs PatchName"))
      mparams->packsName = 
        pmshparams["createPatch Parameters"]["Packs PatchName"].as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Surrounding PatchType"))
      mparams->srrndngPatchType = 
        pmshparams["createPatch Parameters"]["Surrounding PatchType"].as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Packs PatchType"))
      mparams->packsPatchType = 
        pmshparams["createPatch Parameters"]["Packs PatchType"].as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("overwrite?"))
      mparams->_overwritecpMsh = 
        pmshparams["createPatch Parameters"]["overwrite?"].as<bool>();
    if (pmshparams["createPatch Parameters"].contains("Packs Path"))
      mparams->pathPacks = 
        pmshparams["createPatch Parameters"]["Packs Path"].as<std::string>();
    if (pmshparams["createPatch Parameters"].contains("Surrounding Path"))
      mparams->pathSurrounding = 
        pmshparams["createPatch Parameters"]["Surrounding Path"].as<std::string>();
  }
  
  if (pmshparams.contains("foamToSurface Parameters"))
  {
    if (pmshparams["foamToSurface Parameters"].contains("Output File Path"))
      mparams->outSurfName = 
        pmshparams["foamToSurface Parameters"]["Output File Path"].as<std::string>();
  }

  if (pmshparams.contains("surfaceSplitByTopology Parameters"))
  {
    if (pmshparams["surfaceSplitByTopology Parameters"].contains("Input File"))
      mparams->surfFile = 
        pmshparams["surfaceSplitByTopology Parameters"]["Input File"].as<std::string>();
    if (pmshparams["surfaceSplitByTopology Parameters"].contains("Output File"))
      mparams->outSurfFile = 
        pmshparams["surfaceSplitByTopology Parameters"]["Output File"].as<std::string>();
  }

  // creating triSurface directory for workflow
  std::string dir_path = "./constant";
  boost::filesystem::path dir(dir_path);
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating triSurface directory for the snappyHexMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  // directory initialization
  dir = "./constant/triSurface";
  try
  {
    boost::filesystem::create_directory(dir);
  }
  catch (boost::filesystem::filesystem_error &e)
  {
    std::cerr << "Problem in creating triSurface directory for the snappyHexMesh" << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  MeshManipulationFoam* objMsh = new MeshManipulationFoam(mparams);
  const char* nameFile = "a";
  

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
  blockMeshGen* objBM = new blockMeshGen(bmparams);
  objBM->createMeshFromSTL(nameFile);

  // snappyHexMesh
  snappymeshGen* objSHM = new snappymeshGen(snappyparams);
  objSHM->createMeshFromSTL(nameFile);

  // splitMeshRegions
  std::pair<int,int> dirStat = objMsh->splitMshRegions(); // outputs the region number 
                                           // skipped during splitting process
  int skippedDir = dirStat.first;
  int totalRegs = dirStat.second - 1;

  // mergeMeshes
  objMsh->mergeMeshes(skippedDir,nDom);

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

  meshBase* fm = new FOAM::foamMesh(readDB);
  fm->read(regNme);
  vtkMesh* vm = new vtkMesh(fm->getDataSet(),ofname1);
  vm->report();
  vm->write();

  // converts surronding mesh
  regNme = "domain0";
  meshBase* fm2 = new FOAM::foamMesh(readDB);
  fm2->read(regNme);
  vtkMesh* vm2 = new vtkMesh(fm2->getDataSet(),ofname2);
  vm2->report();
  vm2->write();

  // cleaning up
  if (vm)
    delete vm;
  if (fm)
    delete fm;
  if (vm2)
    delete vm2;
  if (fm2)
    delete fm2;
  if (objMsh)
    delete objMsh;
  if (objSHM)
    delete objSHM;
  //if (objCF)
  //  delete objCF;
  if (objBM)
    delete objBM;

  return 0;
}


// TEST macros 
TEST(PackMeshing, Generation)
{
  EXPECT_EQ(0, generate(inp_json));
}

TEST(PackMeshing, NumberOfNodesnCellsPacks)
{
  if (ref)
      delete ref;
  ref = meshBase::Create( inputjson["Pack Reference File"].as<std::string>() );
  meshBase* cmp1 = meshBase::Create( "geom_pack_mesh.vtu" );
  EXPECT_TRUE((cmp1->getNumberOfPoints() == ref->getNumberOfPoints()) &&
              (cmp1->getNumberOfCells() == ref->getNumberOfCells()));
}

TEST(PackMeshing, NumberOfNodesNodesnCellsSurrounding)
{
  if (ref)
      delete ref;
  ref = meshBase::Create( inputjson["Surrounding Reference File"].as<std::string>() );
  meshBase* cmp2 = meshBase::Create( "geom_surrounding_mesh.vtu" );
  EXPECT_TRUE((cmp2->getNumberOfPoints() == ref->getNumberOfPoints()) &&
              (cmp2->getNumberOfCells() == ref->getNumberOfCells()));
}

// test constructor
int main(int argc, char** argv) 
{
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 1);
  inp_json = argv[1];

  if (!inp_json)
  {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }

  // running tests
  int res = RUN_ALL_TESTS();
  
  // clean up
  if (mesh)
      delete mesh;

  return res;
}
