#include <vtkGenericDataObjectReader.h>
#include <vtkSmartPointer.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include "Drivers/MeshGen/CFMeshMeshGenDriver.H"
#include "Drivers/MeshQualityDriver.H"
#include "cfmeshQualityParams.H"
#include "foamMesh.H"
#include "gtest.h"
#include "vtkMesh.H"

const char* inp_json;
const char* inp_json_2;
meshBase* mesh;
FOAM::foamMesh* fmesh;
meshBase* ref;
jsoncons::json inputjson;
jsoncons::json inputjson_2;

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

  auto driver = NEM::DRV::NemDriver::readJSON(inputjson);
  EXPECT_NE(dynamic_cast<NEM::DRV::CFMeshMeshGenDriver *>(driver.get()),
            nullptr);
  driver->execute();

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

  jsoncons::json inputjson_tmp;
  inputStream >> inputjson_tmp;
  inputjson_2 = inputjson_tmp[0];

  auto driver = NEM::DRV::NemDriver::readJSON(inputjson_2);
  EXPECT_NE(dynamic_cast<NEM::DRV::OptimizeMeshQualDriver *>(driver.get()),
            nullptr);
  driver->execute();

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
