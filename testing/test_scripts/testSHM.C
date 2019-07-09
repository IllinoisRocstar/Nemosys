#include "NemDriver.H"
#include "snappymeshGen.H"
#include "snappymeshParams.H"
#include "vtkMesh.H"
#include <gtest.h>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
namespace fs = boost::filesystem;

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

void copyDirectoryRecursively(const fs::path& sourceDir, const fs::path& destinationDir)
{
    if (!fs::exists(sourceDir) || !fs::is_directory(sourceDir))
    {
        throw std::runtime_error("Source directory " + sourceDir.string() + " does not exist or is not a directory");
    }
    if (fs::exists(destinationDir))
    {
        throw std::runtime_error("Destination directory " + destinationDir.string() + " already exists");
    }
    if (!fs::create_directory(destinationDir))
    {
        throw std::runtime_error("Cannot create destination directory " + destinationDir.string());
    }

    for (const auto& dirEnt : fs::recursive_directory_iterator{sourceDir})
    {
        const auto& path = dirEnt.path();
        auto relativePathStr = path.string();
        boost::algorithm::replace_first(relativePathStr, sourceDir.string(), "");
        fs::copy(path, destinationDir / relativePathStr);
    }
}

void snappyLogistics()
{
    const char dir_path[] = "./constant/polyMesh";
    const char cpydir_path[] = "./constant/polyMesh.Orig";
    boost::filesystem::path dir1(dir_path);
	boost::filesystem::path dir2(cpydir_path);
	boost::filesystem::remove_all(dir1);
	copyDirectoryRecursively(dir2,dir1);
	
	return;
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

  snappymeshParams* params = new snappymeshParams();
  std::string ifname = inputjson["Mesh File Options"]
                                ["Input Geometry File"].as<std::string>();  
  std::string ofname = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();  
  jsoncons::json shmparams = inputjson["Meshing Parameters"]["snappyHexMesh Parameters"];

  // required params here
  // cad file
  if (inputjson["Mesh File Options"].has_key("Input Geometry File"))
      params->geomFileName = 
          inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else
  {
      std::cerr << "A geometry file should be supplied.\n";
      throw;        
  }  

        if (shmparams.has_key("SurfPatchName"))
            params->surfRefPatch = shmparams["SurfPatchName"].as<std::string>();
        else
        {
            std::cerr << "A patch name for surface refinement must be defined.\n";
            throw;        
        }

        if (shmparams.has_key("InputPatchName"))
            params->geomPatchName = shmparams["InputPatchName"].as<std::string>();
        else
        {
            std::cerr << "A patch name for input geometry must be defined.\n";
            throw;        
        }

        if (shmparams.has_key("Castellated Mesh"))
          params->_withCastMesh =
            shmparams["Castellated Mesh"].as<bool>();
        if (shmparams.has_key("Snapping"))
          params->_withSnap =
            shmparams["Snapping"].as<bool>();
        if (shmparams.has_key("Layer Addition"))
          params->_withLayers =
            shmparams["Layer Addition"].as<bool>();
        if (shmparams.has_key("CellZones"))
          params->_withCellZones =
            shmparams["CellZones"].as<bool>();
        if (shmparams.has_key("RegionRefine"))
          params->_withGeomRefReg =
            shmparams["RegionRefine"].as<bool>();
        if (shmparams.has_key("SurfaceRefine"))
          params->_withSurfRefReg =
            shmparams["SurfaceRefine"].as<bool>();

        if (shmparams.has_key("maxLocalCells"))
          params->maxLCells =
            shmparams["maxLocalCells"].as<int>();
        if (shmparams.has_key("maxGlobalCells"))
          params->maxGCells =
            shmparams["maxGlobalCells"].as<int>();
        if (shmparams.has_key("minRefCells"))
          params->minRefCells =
            shmparams["minRefCells"].as<int>();
        if (shmparams.has_key("nCellsBetweenLevels"))
          params->cellsBetnLvls =
            shmparams["nCellsBetweenLevels"].as<int>();
        if (shmparams.has_key("surfaceRefinementLvlMin"))
          params->refSurfLvlMin =
            shmparams["surfaceRefinementLvlMin"].as<int>();
        if (shmparams.has_key("surfaceRefinementLvlMax"))
          params->refSurfLvlMax =
            shmparams["surfaceRefinementLvlMax"].as<int>();
        if (shmparams.has_key("resolveFeatureAngle"))
          params->featAngle =
            shmparams["resolveFeatureAngle"].as<double>();
        if (shmparams.has_key("locationInMeshX"))
          params->locMeshX =
            shmparams["locationInMeshX"].as<double>();
        if (shmparams.has_key("locationInMeshY"))
          params->locMeshY =
            shmparams["locationInMeshY"].as<double>();
        if (shmparams.has_key("locationInMeshZ"))
          params->locMeshZ =
            shmparams["locationInMeshZ"].as<double>();
        if (shmparams.has_key("allowFreeStandingZoneFaces"))
          params->_alwFreeZone =
            shmparams["allowFreeStandingZoneFaces"].as<double>();
        if (shmparams.has_key("nSmoothPatch"))
          params->snapSmthPatch =
            shmparams["nSmoothPatch"].as<int>();
        if (shmparams.has_key("tolerance"))
          params->snapTol =
            shmparams["tolerance"].as<double>();
        if (shmparams.has_key("snapSolveIter"))
          params->solveSnapIter =
            shmparams["snapSolveIter"].as<int>();
        if (shmparams.has_key("snapRelaxIter"))
          params->relaxSnapIter =
            shmparams["snapRelaxIter"].as<int>();
        if (shmparams.has_key("relativeSizes"))
          params->_relSize =
            shmparams["relativeSizes"].as<double>();
        if (shmparams.has_key("expansionRatio"))
          params->expRatio =
            shmparams["expansionRatio"].as<double>();
        if (shmparams.has_key("finalLayerThickness"))
          params->finLThick =
            shmparams["finalLayerThickness"].as<double>();
        if (shmparams.has_key("minThickness"))
          params->minThick =
            shmparams["minThickness"].as<double>();
        if (shmparams.has_key("nGrow"))
          params->nGrow =
            shmparams["nGrow"].as<int>();
        if (shmparams.has_key("featureAngle"))
          params->lyrFeatAngle =
            shmparams["featureAngle"].as<double>();
        if (shmparams.has_key("nRelaxIter"))
          params->lyrRelaxIter =
            shmparams["nRelaxIter"].as<int>();
        if (shmparams.has_key("nSmoothSurfaceNormals"))
          params->lyrSmthSurfNorm =
            shmparams["nSmoothSurfaceNormals"].as<int>();
        if (shmparams.has_key("nSmoothNormals"))
          params->lyrSmthNorm =
            shmparams["nSmoothNormals"].as<int>();
        if (shmparams.has_key("nSmoothThickness"))
          params->lyrSmthThick =
            shmparams["nSmoothThickness"].as<int>();
        if (shmparams.has_key("maxFaceThicknessRatio"))
          params->lyrMaxFcTR =
            shmparams["maxFaceThicknessRatio"].as<double>();
        if (shmparams.has_key("maxThicknessToMedialRatio"))
          params->lyrMaxThickTMR =
            shmparams["maxThicknessToMedialRatio"].as<double>();
        if (shmparams.has_key("minMedialAxisAngle"))
          params->lyrMinMedAngl =
            shmparams["minMedialAxisAngle"].as<double>();
        if (shmparams.has_key("nBufferCellsNoExtrude"))
          params->lyrBuffrCells =
            shmparams["nBufferCellsNoExtrude"].as<int>();
        if (shmparams.has_key("nLayerIter"))
          params->lyrIter =
            shmparams["nLayerIter"].as<int>();
        if (shmparams.has_key("maxNonOrtho"))
          params->qcMaxNOrtho =
            shmparams["maxNonOrtho"].as<int>();
        if (shmparams.has_key("maxBoundarySkewness"))
          params->qcMaxBndrySkew =
            shmparams["maxBoundarySkewness"].as<double>();
        if (shmparams.has_key("maxInternalSkewness"))
          params->qcMaxIntSkew =
            shmparams["maxInternalSkewness"].as<double>();
        if (shmparams.has_key("maxConcave"))
          params->qcMaxConc =
            shmparams["maxConcave"].as<double>();
        if (shmparams.has_key("minVol"))
          params->qcMinVol =
            shmparams["minVol"].as<double>();
        if (shmparams.has_key("minTetQuality"))
          params->qcMinTetQ =
            shmparams["minTetQuality"].as<double>();
        if (shmparams.has_key("minArea"))
          params->qcMinArea =
            shmparams["minArea"].as<double>();
        if (shmparams.has_key("minTwist"))
          params->qcMinTwist =
            shmparams["minTwist"].as<double>();
        if (shmparams.has_key("minFaceWeight"))
          params->qcMinFaceW =
            shmparams["minFaceWeight"].as<double>();
        if (shmparams.has_key("minVolRatio"))
          params->qcMinVolRto =
            shmparams["minVolRatio"].as<double>();
        if (shmparams.has_key("minDeterminant"))
          params->qcMinDet =
            shmparams["minDeterminant"].as<double>();
        if (shmparams.has_key("minTriangleTwist"))
          params->qcMinTrTwist =
            shmparams["minTriangleTwist"].as<double>();
        if (shmparams.has_key("qcnSmoothScale"))
          params->qcSmthScale =
            shmparams["qcnSmoothScale"].as<int>();
        if (shmparams.has_key("errorReduction"))
          params->qcErrRedctn =
            shmparams["errorReduction"].as<double>();
        if (shmparams.has_key("mergeTolerance"))
          params->mergeTol =
            shmparams["mergeTolerance"].as<double>();


        std::string cap2 = "GeomRefinementRegions";
        if (shmparams.has_key(cap2))
        {
          params->_withGeomRefReg = 1;

          for (auto jptch2 : shmparams[cap2].array_range())
          {
            shmGeomRefine geomRef;

            if (jptch2.has_key("PatchName"))
              geomRef.patchNm = jptch2["PatchName"].as<std::string>();

            if (jptch2.has_key("searchableShape"))
              geomRef.searchableName = jptch2["searchableShape"].as<std::string>();

            if (jptch2.has_key("shapeParams1"))
              geomRef.shapeParameters1 = jptch2["shapeParams1"].as<std::string>();

            if (jptch2.has_key("shapeParams2"))
              geomRef.shapeParameters2 = jptch2["shapeParams2"].as<std::string>();

            if (jptch2.has_key("Radius"))
              geomRef.rad = jptch2["Radius"].as<double>();

            if (jptch2.has_key("Mode"))
              geomRef.mode = jptch2["Mode"].as<std::string>();

            if (jptch2.has_key("MinLevel"))
              geomRef.minLvl = jptch2["MinLevel"].as<int>();

            if (jptch2.has_key("MaxLevel"))
              geomRef.maxLvl = jptch2["MaxLevel"].as<int>();

            (params->geomRefs).push_back(geomRef);

          }
        }


        std::string cap3 = "SurfaceRefinementRegions";
        if (shmparams.has_key(cap3))
        {
          params->_withSurfRefReg = 1;

          for (auto jptch3 : shmparams[cap3].array_range())
          {

            shmRegionRef surfRef;

            if (jptch3.has_key("PatchName"))
              surfRef.refPatchNm = jptch3["PatchName"].as<std::string>();

            if (jptch3.has_key("MinLevel"))
              surfRef.minLvl = jptch3["MinLevel"].as<int>();

            if (jptch3.has_key("MaxLevel"))
              surfRef.maxLvl = jptch3["MaxLevel"].as<int>();

            (params->surfRefs).push_back(surfRef);
          }
        }

  snappymeshGen* generator = new snappymeshGen(dynamic_cast<snappymeshParams*>(params));
  generator->createMeshFromSTL("");
  mesh = vtkMesh::Create(generator->getDataSet(), ofname);
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();

  return 0;
}

// Write a new name for snappyHexMesh
// TEST macros 
TEST(snappyHexMesh, Generation)
{
  EXPECT_EQ(0, generate(inp_json));
}

TEST(snappyHexMesh, NumberOfNodes)
{
  if (ref)
      delete ref;
  ref = meshBase::Create( inputjson["Reference File"].as<std::string>() );
  EXPECT_EQ( mesh->getNumberOfPoints(), ref->getNumberOfPoints() );
}

TEST(snappyHexMesh, NumberOfCells)
{
  if (ref)
      delete ref;
  ref = meshBase::Create( inputjson["Reference File"].as<std::string>() );
  EXPECT_EQ( mesh->getNumberOfCells(), ref->getNumberOfCells() );
}

// test constructor
int main(int argc, char** argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 1);
  inp_json = argv[1];

  if (!inp_json)
  {
  	std::cerr << "No input file defined" << std::endl;
  }

  int res = RUN_ALL_TESTS();
  
  // clean up
  if (mesh)
      delete mesh;
  
  snappyLogistics(); // clean up and replace mesh

  return res;


}
