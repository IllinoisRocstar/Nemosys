#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <string>
#include <Drivers/NemDriver.H>
#include <Drivers/PackMesh/PackMeshDriver.H>
#include <Mesh/geoMeshFactory.H>
#include <Mesh/meshBase.H>

const char *inp_json;
NEM::MSH::geoMeshBase *mesh;
NEM::MSH::geoMeshBase *ref;
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

  // Call packmesh readjson
  auto pckmshdrvObj = NEM::DRV::NemDriver::readJSON(inputjson);
  pckmshdrvObj->execute();

  return 0;
}

// TEST macros
TEST(PackMeshing, Generation) { EXPECT_EQ(0, generate(inp_json)); }

TEST(PackMeshing, NumberOfNodesNodesnCellsMerged) {
  if (ref) delete ref;
  ref = NEM::MSH::Read(inputjson["Merged Reference File"].as<std::string>());
  auto *cmp2 = NEM::MSH::Read("packmesh.vtu");
  EXPECT_TRUE(std::fabs(cmp2->getNumberOfPoints() / ref->getNumberOfPoints() -
                        1) < 0.03);
  EXPECT_TRUE(
      std::fabs(cmp2->getNumberOfCells() / ref->getNumberOfCells() - 1) < 0.03);
  cmp2->Delete();
  ref->Delete();
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
