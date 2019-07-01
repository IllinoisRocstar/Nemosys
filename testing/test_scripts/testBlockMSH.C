#include "NemDriver.H"
#include "blockMeshGen.H"
#include "blockMeshParams.H"
#include "vtkMesh.H"
#include <gtest.h>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

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

  blockMeshParams* params = new blockMeshParams();
  std::string ifname = inputjson["Mesh File Options"]
                                ["Input Geometry File"].as<std::string>();  
  std::string ofname = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();  
  jsoncons::json bmshparams = 
        inputjson["Meshing Parameters"]["blockMesh Parameters"];

  // required params here
  // cad file
  if (inputjson["Mesh File Options"].has_key("Input Dict File"))
      params->_ownBlockMshDict = 
        inputjson["Mesh File Options"]["Input Dict File"].as<bool>();

    
  // Parameter parsing starts here
  if (bmshparams.has_key("Block Geometry"))
    params->_isBlock =
      bmshparams["Block Geometry"].as<bool>();
  if (bmshparams.has_key("Sphere Geometry"))
    params->_isSphere =
      bmshparams["Sphere Geometry"].as<bool>();
  if (bmshparams.has_key("Cylinder/Tapered_Cone Geometry"))
    params->_isCylinder_TCone =
      bmshparams["Cylinder/Tapered_Cone Geometry"].as<bool>();
  if (bmshparams.has_key("scaleToMeters"))
    params->cnvrtToMeters = 
    bmshparams["scaleToMeters"].as<double>();
  if (bmshparams.has_key("XdirectionCells"))
    params->cellsXDir = 
    bmshparams["XdirectionCells"].as<int>();
  if (bmshparams.has_key("YdirectionCells"))
    params->cellsYDir = 
    bmshparams["YdirectionCells"].as<int>();
  if (bmshparams.has_key("ZdirectionCells"))
    params->cellsZDir = 
    bmshparams["ZdirectionCells"].as<int>();

  if (bmshparams.has_key("Cell_Size")){
    params->_cellSizeDefined = true;
    params->cellSize = 
     bmshparams["Cell_Size"].as<double>();
  }
  else
  {
    params->cellSize = -1;
    params->_cellSizeDefined = false;
  }
      
  if (bmshparams.has_key("Block Parameters"))
  {

    if (bmshparams["Block Parameters"].has_key("Auto_Generate"))
    {
      params->_autoGenerateBox = true;

      if (inputjson["Mesh File Options"].has_key("Input Geometry File"))
        params->packFileName = 
            inputjson["Mesh File Options"]
                    ["Input Geometry File"].as<std::string>();
      else
      {
        std::cerr << "A geometry file should be supplied.\n";
        throw;        
      }

      if (bmshparams["Block Parameters"]
          ["Auto_Generate"].has_key("Offset_XDir"))
            params->offsetX = 
            bmshparams["Block Parameters"]
              ["Auto_Generate"]["Offset_XDir"].as<double>();
      else{
        params->offsetX = 0.1;
      }
      if (bmshparams["Block Parameters"]
          ["Auto_Generate"].has_key("Offset_YDir"))
           params->offsetY = 
            bmshparams["Block Parameters"]
            ["Auto_Generate"]["Offset_YDir"].as<double>();
      else{
        params->offsetY = 0.1;
      }
      if (bmshparams["Block Parameters"]
          ["Auto_Generate"].has_key("Offset_ZDir"))
           params->offsetZ = 
            bmshparams["Block Parameters"]
            ["Auto_Generate"]["Offset_ZDir"].as<double>();
      else{
        params->offsetZ = 0.1;
      }

    }
    else{
      params->_autoGenerateBox = false;
    }
      
    if (bmshparams["Block Parameters"].has_key("X1"))
      params->initX = bmshparams["Block Parameters"]["X1"].as<double>();
    if (bmshparams["Block Parameters"].has_key("Y1"))
      params->initY = bmshparams["Block Parameters"]["Y1"].as<double>();
    if (bmshparams["Block Parameters"].has_key("Z1"))
      params->initZ = bmshparams["Block Parameters"]["Z1"].as<double>();
    if (bmshparams["Block Parameters"].has_key("LengthX"))
      params->lenX = bmshparams["Block Parameters"]["LengthX"].as<double>();
    if (bmshparams["Block Parameters"].has_key("LengthY"))
      params->lenY = bmshparams["Block Parameters"]["LengthY"].as<double>();
    if (bmshparams["Block Parameters"].has_key("LengthZ"))
      params->lenZ = bmshparams["Block Parameters"]["LengthZ"].as<double>();
    if (bmshparams["Block Parameters"].has_key("GradingXdir"))
      params->smplGradingX = 
        bmshparams["Block Parameters"]["GradingXdir"].as<double>();
    if (bmshparams["Block Parameters"].has_key("GradingYdir"))
      params->smplGradingY = 
        bmshparams["Block Parameters"]["GradingYdir"].as<double>();
    if (bmshparams["Block Parameters"].has_key("GradingZdir"))
      params->smplGradingZ = 
        bmshparams["Block Parameters"]["GradingZdir"].as<double>();
  }
    
  if (bmshparams.has_key("Sphere Parameters"))
  {
      
    if (bmshparams["Sphere Parameters"].has_key("Center Y"))
      params->centerX =
        bmshparams["Sphere Parameters"]["Center X"].as<double>();
    if (bmshparams["Sphere Parameters"].has_key("Center Y"))
      params->centerY =
        bmshparams["Sphere Parameters"]["Center Y"].as<double>();
    if (bmshparams["Sphere Parameters"].has_key("Center Z"))
      params->centerZ =
        bmshparams["Sphere Parameters"]["Center Z"].as<double>();
    if (bmshparams["Sphere Parameters"].has_key("Radius"))
      params->radius =
        bmshparams["Sphere Parameters"]["Radius"].as<double>();
    if (bmshparams["Sphere Parameters"].has_key("GradingXdir"))
      params->sphrGradingX =
        bmshparams["Sphere Parameters"]["GradingXdir"].as<double>();
    if (bmshparams["Sphere Parameters"].has_key("GradingYdir"))
      params->sphrGradingY =
        bmshparams["Sphere Parameters"]["GradingYdir"].as<double>();
    if (bmshparams["Sphere Parameters"].has_key("GradingXdir"))
      params->sphrGradingZ =
        bmshparams["Sphere Parameters"]["GradingZdir"].as<double>();
  }
  
  if (bmshparams.has_key("Cylinder/Tapered_Cone Parameters"))
  {
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Center X"))
      params->centerCyl[0] =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["Center X"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Center Y"))
      params->centerCyl[1] =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Y"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Center Z"))
      params->centerCyl[2] =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["Center Z"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Radius1"))
      params->radius1 =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"].as<double>();
      
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Radius2")){
      params->radius2 =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius2"].as<double>();
    }
    else{
      params->radius2 =
          bmshparams["Cylinder/Tapered_Cone Parameters"]["Radius1"].as<double>();
    }
      
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("GradingXdir"))
      params->cylGrading[0] =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingXdir"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("GradingYdir"))
      params->cylGrading[1] =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingYdir"].as<double>();
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("GradingXdir"))
      params->cylGrading[2] =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["GradingZdir"].as<double>();
      
    if (bmshparams["Cylinder/Tapered_Cone Parameters"].has_key("Height"))
      params->height =
        bmshparams["Cylinder/Tapered_Cone Parameters"]["Height"].as<double>();
  }

  auto *generator =
        new blockMeshGen(dynamic_cast<blockMeshParams *>(params));
  generator->createMeshFromSTL("");
  mesh = vtkMesh::Create(generator->getDataSet(), ofname);
  mesh->setFileName(ofname);
  mesh->report();
  mesh->write();

  return 0;
}

// TEST macros 
TEST(blockMesh, Generation)
{
  EXPECT_EQ(0, generate(inp_json));
}

TEST(blockMesh, NumberOfNodes)
{
  if (ref)
      delete ref;
  ref = meshBase::Create( inputjson["Reference File"].as<std::string>() );
  EXPECT_EQ( mesh->getNumberOfPoints(), ref->getNumberOfPoints() );
}

TEST(blockMesh, NumberOfCells)
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
    throw;
  }

  int res = RUN_ALL_TESTS();
  
  // clean up
  if (mesh)
      delete mesh;

  return res;


}
