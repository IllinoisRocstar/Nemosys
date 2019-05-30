#include "NemDriver.H"
#include "cfmeshGen.H"
#include "cfmeshParams.H"
#include "cfmeshGen.H"
#include "cfmeshQualityParams.H"
#include "MeshQuality.H"
#include "foamMesh.H"
#include "vtkMesh.H"
#include "gtest.h"
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <vtkGenericDataObjectReader.h>
#include <vtkSmartPointer.h>

const char* inp_json;
const char* inp_json_2;
meshBase* mesh;
FOAM::foamMesh* fmesh;
meshBase* ref;
json inputjson;
json inputjson_2;

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
  
  json inputjson_tmp;
  inputStream >> inputjson_tmp;
  inputjson = inputjson_tmp[0];  

  cfmeshParams* params = new cfmeshParams();
  std::string ifname = inputjson["Mesh File Options"]
                                ["Input Geometry File"].as<std::string>();  
  std::string ofname = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();  
  json cfmparams = inputjson["Meshing Parameters"]["CFMesh Parameters"];

  // required params here
  // cad file
  if (inputjson["Mesh File Options"].has_key("Input Geometry File"))
      params->geomFilePath = 
          inputjson["Mesh File Options"]["Input Geometry File"].as<std::string>();
  else
  {
      std::cerr << "A geometry file should be supplied.\n";
      throw;        
  }
  
  // mesh generator
  if (cfmparams.has_key("Generator"))
      params->generator = cfmparams["Generator"].as<std::string>();
  else
  {
      std::cerr << "A mesh generation method should be selected.\n";
      std::cerr << "Options: cartesian2D tetMesh\n";
      throw;        
  }
  
  // rest of params are optional
  if (cfmparams.has_key("MaxCellSize"))
    params->maxCellSize = 
        cfmparams["MaxCellSize"].as<double>();
  if (cfmparams.has_key("MinCellSize"))
    params->minCellSize = 
        cfmparams["MinCellSize"].as<double>();
  if (cfmparams.has_key("BoundaryCellSize"))
    params->bndryCellSize = 
        cfmparams["BoundaryCellSize"].as<double>();
  if (cfmparams.has_key("KeepCellsIntersectingBoundary"))
    params->keepCellIB = 
        cfmparams["KeepCellsIntersectingBoundary"].as<double>();
  if (cfmparams.has_key("CheckForGluedMesh"))
    params->chkGluMsh = 
        cfmparams["CheckForGluedMesh"].as<double>();

  // optional capability
  std::string cap = "BoundaryLayers";
  if (cfmparams.has_key(cap))
  {
      params->_withBndLyr = true;
      params->blNLyr = cfmparams[cap]["NLayers"].as<double>();
      params->blThkRto = cfmparams[cap]["ThicknessRatio"].as<double>();
      if (cfmparams[cap].has_key("MaxFirstLayerThickness"))
          params->maxFrstLyrThk = cfmparams[cap]["MaxFirstLayerThickness"].as<double>();
      if (cfmparams[cap].has_key("AllowDiscontinuity"))
          params->alwDiscont = cfmparams[cap]["AllowDiscontinuity"].as<bool>();

      // patch boundary layers
      std::string subcap = "PatchBoundaryLayers";
      if (cfmparams[cap].has_key(subcap))
      {
          params->_withBndLyrPtch = true;
          for (auto jptch : cfmparams[cap][subcap].array_range())
          {
              cfmPtchBndLyr blPatch;
              blPatch.patchName = jptch["PatchName"].as<std::string>();
              if (jptch.has_key("AllowDiscontinuity"))
                  blPatch.alwDiscont = jptch["AllowDiscontinuity"].as<bool>();
              else
                  blPatch.alwDiscont = false;
              if (jptch.has_key("MaxFirstLayerThickness"))
                  blPatch.maxFrstLyrThk = jptch["MaxFirstLayerThickness"].as<int>();
              else
                  blPatch.maxFrstLyrThk = -1;
              if (jptch.has_key("NLayers"))
                  blPatch.blNLyr = jptch["NLayers"].as<int>();
              else
                  blPatch.blNLyr = -1;
              if (jptch.has_key("ThicknessRatio"))
                  blPatch.blThkRto = jptch["ThicknessRatio"].as<double>();
              else
                  blPatch.blThkRto = -1.;
              (params->blPatches).push_back(blPatch);
          }
                      
      }
  }

  // optional capability
  cap = "SurfaceFeatureEdges";
  if (cfmparams.has_key(cap))
  {
      params->_withSrfEdg = true;
      params->srfEdgAng = cfmparams[cap]["Angle"].as<double>();
  }

  // optional capability
  cap = "ObjectRefinements";
  if (cfmparams.has_key(cap))
  {
      params->_withObjRfn = true;
      for (auto refObj : cfmparams[cap].array_range())
      {
          cfmObjRef objRef;
          objRef.name= refObj["Name"].as<std::string>();
          for (const auto&  prm: refObj["Params"].object_range())
          {
              std::string key = prm.key();
              std::string val = prm.value().as<std::string>();
              objRef.params[key] = val;
          }
          (params->objRefLst).push_back(objRef);
      }
  }

  // optional capability
  cap = "ImproveMeshQuality";
  if (cfmparams.has_key(cap))
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
  if (cfmparams.has_key(cap))
  {
      params->_withLclRef = true;
      for (auto jptch : cfmparams[cap].array_range())
      {
          cfmLclRefPatch refPatch;
          refPatch.patchName = jptch["PatchName"].as<std::string>();
          if (jptch.has_key("AdditionalRefinementLevels"))
              refPatch.aditRefLvls = jptch["AdditionalRefinementLevels"].as<int>();
          else
              refPatch.aditRefLvls = -1;
          if (jptch.has_key("CellSize"))
              refPatch.cellSize = jptch["CellSize"].as<double>();
          else
              refPatch.cellSize = -1.;
          (params->refPatches).push_back(refPatch);
      }
  }

  // optional capability
  cap = "RenameBoundary";
  if (cfmparams.has_key(cap))
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
  delete mesh;

  return 0;
}



int optimize(const char* jsonF)
{
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if(!inputStream.good())
  {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }
  
  json inputjson_tmp;
  inputStream >> inputjson_tmp;
  inputjson_2 = inputjson_tmp[0];
  
  cfmshQualityParams* params = new cfmshQualityParams();
  MeshQuality* mq = new MeshQuality(params);
  
  // carry out schedules
  for (auto jsch : inputjson_2["Schedule"].array_range())
  {
      std::string qImpMtd = jsch["Method"].as<std::string>(); 
      if (!qImpMtd.compare("meshOptimizer"))
      {
        params->nIterations = 
            jsch["Params"]["NIterations"].as<int>();
        params->nLoops = jsch["Params"]["NLoops"].as<int>();
        params->qualThrsh = 
            jsch["Params"]["QualityThreshold"].as<double>();
        params->nSrfItr = 
            jsch["Params"]["NSurfaceIterations"].as<int>();
        if (jsch["Params"].has_key("ConstrainedCellSet"))
        {
            params->_withConstraint = true;
            params->consCellSet = 
                jsch["Params"]["ConstrainedCellSet"].as<std::string>();
        }

        mq->cfmOptimize();
      }
  }

  fmesh = new FOAM::foamMesh(true);

  return 0;
}


// TEST macros
TEST(CfMesh, Generation)
{
  EXPECT_EQ(0, generate(inp_json));
}

TEST(CfMesh, Optimization)
{
  EXPECT_EQ(0, optimize(inp_json_2));
}

TEST(CfMesh, NumberOfNodes)
{
  // reading legacy vtk files
  vtkSmartPointer<vtkGenericDataObjectReader> reader =
      vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader->SetFileName( 
          inputjson["Reference File"].as<std::string>().c_str());
  reader->Update();

  EXPECT_EQ( fmesh->getNumberOfPoints(), 
          reader->GetUnstructuredGridOutput()->GetNumberOfPoints() );
}

TEST(CfMesh, NumberOfCells)
{
  // reading legacy vtk files
  vtkSmartPointer<vtkGenericDataObjectReader> reader =
      vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader->SetFileName( 
          inputjson["Reference File"].as<std::string>().c_str());
  reader->Update();

  EXPECT_EQ( fmesh->getNumberOfCells(), 
          reader->GetUnstructuredGridOutput()->GetNumberOfCells() );
}

// test constructor
int main(int argc, char** argv) 
{
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc >= 3);
  inp_json = argv[1];
  inp_json_2 = argv[2];
  
  // run tests
  int res = RUN_ALL_TESTS();
  
  // clean up
  //if (mesh)
  //    delete mesh;
  //if (fmesh)
  //    delete fmesh;

  return res;
}
