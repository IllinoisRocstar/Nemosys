#include "NemDriver.H"
#include "cfmeshGen.H"
#include "cfmeshParams.H"
#include "vtkMesh.H"
#include <gtest.h>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>

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

  cfmeshParams* params = new cfmeshParams();
  std::string ifname = inputjson["Mesh File Options"]
                                ["Input Geometry File"].as<std::string>();  
  std::string ofname = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();
  jsoncons::json cfmparams = inputjson["Meshing Parameters"]["CFMesh Parameters"];

  // required params here
  // cad file
  if (inputjson["Mesh File Options"].contains("Input Geometry File"))
      params->geomFilePath = 
          inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else
  {
      std::cerr << "A geometry file should be supplied.\n";
      throw;        
  }
  
  // mesh generator
  if (cfmparams.contains("Generator"))
      params->generator = cfmparams["Generator"].as<std::string>();
  else
  {
      std::cerr << "A mesh generation method should be selected.\n";
      std::cerr << "Options: cartesian2D tetMesh\n";
      throw;        
  }
  
  // rest of params are optional
  if (cfmparams.contains("MaxCellSize"))
    params->maxCellSize = 
        cfmparams["MaxCellSize"].as<double>();
  if (cfmparams.contains("MinCellSize"))
    params->minCellSize = 
        cfmparams["MinCellSize"].as<double>();
  if (cfmparams.contains("BoundaryCellSize"))
    params->bndryCellSize = 
        cfmparams["BoundaryCellSize"].as<double>();
  if (cfmparams.contains("KeepCellsIntersectingBoundary"))
    params->keepCellIB = 
        cfmparams["KeepCellsIntersectingBoundary"].as<double>();
  if (cfmparams.contains("CheckForGluedMesh"))
    params->chkGluMsh = 
        cfmparams["CheckForGluedMesh"].as<double>();

  // optional capability
  std::string cap = "BoundaryLayers";
  if (cfmparams.contains(cap))
  {
      params->_withBndLyr = true;
      params->blNLyr = cfmparams[cap]["NLayers"].as<double>();
      params->blThkRto = cfmparams[cap]["ThicknessRatio"].as<double>();
      if (cfmparams[cap].contains("MaxFirstLayerThickness"))
          params->maxFrstLyrThk = cfmparams[cap]["MaxFirstLayerThickness"].as<double>();
      if (cfmparams[cap].contains("AllowDiscontinuity"))
          params->alwDiscont = cfmparams[cap]["AllowDiscontinuity"].as<bool>();

      // patch boundary layers
      std::string subcap = "PatchBoundaryLayers";
      if (cfmparams[cap].contains(subcap))
      {
          params->_withBndLyrPtch = true;
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
              (params->blPatches).push_back(blPatch);
          }
                      
      }
  }

  // optional capability
  cap = "SurfaceFeatureEdges";
  if (cfmparams.contains(cap))
  {
      params->_withSrfEdg = true;
      params->srfEdgAng = cfmparams[cap]["Angle"].as<double>();
  }

  // optional capability
  cap = "ObjectRefinements";
  if (cfmparams.contains(cap))
  {
      params->_withObjRfn = true;
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
          (params->objRefLst).push_back(objRef);
      }
  }

  // optional capability
  cap = "ImproveMeshQuality";
  if (cfmparams.contains(cap))
  {
      params->_withMshQlt = true;
      params->qltNItr = cfmparams[cap]["NIterations"].as<int>();
      params->qltNLop = cfmparams[cap]["NLoops"].as<int>();
      params->qltQltThr = cfmparams[cap]["QualityThreshold"].as<double>();
      params->qltNSrfItr = cfmparams[cap]["NSurfaceIterations"].as<int>();
      params->qltConCelSet = 
          cfmparams[cap].get_with_default("ConstrainedCellsSet","none");
  }

  // optional capability
  cap = "LocalRefinement";
  if (cfmparams.contains(cap))
  {
      params->_withLclRef = true;
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
          (params->refPatches).push_back(refPatch);
      }
  }

  // optional capability
  cap = "RenameBoundary";
  if (cfmparams.contains(cap))
  {
      params->_withRenBndry = true;
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
      (params->renBndry) = renBndry;
  }

  cfmeshGen* generator = new cfmeshGen(dynamic_cast<cfmeshParams*>(params));
  generator->createMeshFromSTL("");
  mesh = vtkMesh::Create(generator->getDataSet(), ofname);
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();

  return 0;
}


// TEST macros
TEST(CfMesh, Generation)
{
  EXPECT_EQ(0, generate(inp_json));
}

TEST(CfMesh, NumberOfNodes)
{
  if (ref)
      delete ref;
  ref = meshBase::Create( inputjson["Reference File"].as<std::string>() );
  EXPECT_EQ( mesh->getNumberOfPoints(), ref->getNumberOfPoints() );
}

TEST(CfMesh, NumberOfCells)
{
  if (ref)
      delete ref;
  ref = meshBase::Create( inputjson["Reference File"].as<std::string>() );
  EXPECT_EQ( mesh->getNumberOfCells(), ref->getNumberOfCells() );
}

// NOTE: This test is unreliable for binary files.
//TEST(CfMesh, FileDiff)
//{
//  bool res = compareFiles
//      ( 
//       inputjson["Reference File"].as<std::string>(),
//       inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>() 
//      );
//  EXPECT_EQ(res, 1);
//}


// test constructor
int main(int argc, char** argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 2);
  inp_json = argv[1];
  
  // run tests
  int res = RUN_ALL_TESTS();
  
  // clean up
  if (mesh)
      delete mesh;

  return res;
}
