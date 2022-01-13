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
#include <gtest/gtest.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <Drivers/MeshGen/CFMeshMeshGenDriver.H>
#include <MeshGeneration/cfmeshGen.H>
#include <MeshGeneration/cfmeshParams.H>
#include <Mesh/geoMeshFactory.H>
#include <Mesh/vtkMesh.H>

const char* inp_json;
NEM::MSH::geoMeshBase* mesh;
NEM::MSH::geoMeshBase* ref;
jsoncons::json inputjson;

// Aux functions
bool compareFiles(const std::string& p1, const std::string& p2) {
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
int generate(const char* jsonF) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson_tmp;
  inputStream >> inputjson_tmp;
  inputjson = inputjson_tmp[0];

  auto driver = std::unique_ptr<NEM::DRV::CFMeshMeshGenDriver>(
      dynamic_cast<NEM::DRV::CFMeshMeshGenDriver*>(
          NEM::DRV::NemDriver::readJSON(inputjson).release()));
  EXPECT_NE(driver, nullptr);

  // Copied from CFMeshMeshGenDriver::execute()
  auto paramsCopy = driver->getParams();
  cfmeshGen generator{&paramsCopy};
  // Parameter not used
  generator.createMeshFromSTL(nullptr);
  const auto& inputGeoFile = driver->getFiles().inputGeoFile;
  std::string newname =
      inputGeoFile.substr(0, inputGeoFile.find_last_of('.')) + ".vtu";

  auto* mshWriter = NEM::MSH::Read(".foam");
  mesh = NEM::MSH::New(driver->getFiles().outputMeshFile);
  mesh->takeGeoMesh(mshWriter);
  mesh->write(driver->getFiles().outputMeshFile);
  mshWriter->Delete();

  return 0;
}

// TEST macros
TEST(CfMesh, Generation) { EXPECT_EQ(0, generate(inp_json)); }

TEST(CfMesh, NumberOfNodes) {
  if (ref) delete ref;
  ref = NEM::MSH::Read(inputjson["Reference File"].as<std::string>());
  EXPECT_EQ(mesh->getNumberOfPoints(), ref->getNumberOfPoints());
}

TEST(CfMesh, NumberOfCells) {
  EXPECT_EQ(mesh->getNumberOfCells(), ref->getNumberOfCells());
}

// NOTE: This test is unreliable for binary files.
// TEST(CfMesh, FileDiff)
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
  mesh->Delete();
  ref->Delete();

  return res;
}
