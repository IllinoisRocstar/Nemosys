#include <fstream>

#include <gtest.h>

#include "NemDriver.H"
#include "meshBase.H"

const char *inp_json;
const char *atr_example_REF;

// Test implementations
std::string nucmeshGen(const char *jsonF) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  if (inputjson.contains("Geometry and Mesh")) {
    std::unique_ptr<NemDriver> nemdrvobj =
        std::unique_ptr<NemDriver>(NemDriver::readJSON(inputjson));
  }

  std::string ifname = "atr_example.msh";

  return ifname;
}

// TEST macros
TEST(NucMesh, AtrExample) {
  int numNodes1, numNodes2, numCells1, numCells2, ret1, ret2;
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(nucmeshGen(inp_json));
  numNodes1 = mesh->getNumberOfPoints();
  numCells1 = mesh->getNumberOfCells();

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(atr_example_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  numCells2 = refMesh->getNumberOfCells();

  if (numNodes1 == numNodes2)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);

  // if (numCells1 == numCells2)
  //   ret2 = 0;
  // else
  //   ret2 = 1;
  // ASSERT_EQ(0, ret2);
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 3);
  inp_json = argv[1];
  atr_example_REF = argv[2];

  if (!inp_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
