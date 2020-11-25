#include <NemDriver.H>
#include <exoMesh.H>
#include <gtest.h>
#include <fstream>

const char *bench1_json;
const char *bench1_ref;
const char *bench2_json;
const char *bench2_ref;
NEM::MSH::EXOMesh::exoMesh *rm;
NEM::MSH::EXOMesh::exoMesh *nm;

// Aux functions
int compareFiles(std::string fn1, std::string fn2) {
  fstream f1, f2;
  int pass = 1;
  f1.open(fn1, fstream::in);
  f2.open(fn2, fstream::in);
  char string1[256], string2[256];
  int j = 0;
  while (!f1.eof()) {
    f1.getline(string1, 256);
    f2.getline(string2, 256);
    j++;
    if (strcmp(string1, string2) != 0) {
      cout << j << "-the strings are not equal"
           << "\n";
      cout << "   " << string1 << "\n";
      cout << "   " << string2 << "\n";
      pass = 0;
      break;
    }
  }
  return pass;
}

// Test implementations
int exoTestConv(const char *jsonF, const char *ofname, const char *refname) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  for (const auto &prog : inputjson.array_range()) {
    std::unique_ptr<NEM::DRV::NemDriver> nemdrvobj =
        std::unique_ptr<NEM::DRV::NemDriver>(
            NEM::DRV::NemDriver::readJSON(prog));
  }

  rm = new NEM::MSH::EXOMesh::exoMesh(std::string(refname));
  nm = new NEM::MSH::EXOMesh::exoMesh(ofname);
  rm->read();
  nm->read();

  return 0;
}

int exoTestNodeNum() {
  return (rm->getNumberOfNodes() == nm->getNumberOfNodes());
}

int exoTestNumElm() {
  return (rm->getNumberOfElements() == nm->getNumberOfElements());
}

int exoTestNodeSetNum() {
  return (rm->getNumberOfNodeSets() == nm->getNumberOfNodeSets());
}

int exoTestElmBlkNum() {
  return (rm->getNumberOfElementBlocks() == nm->getNumberOfElementBlocks());
}

int exoTestNodeSetNames() {
  int nns = rm->getNumberOfNodeSets();
  for (int ins = 1; ins < nns; ins++)
    if (rm->getNdeSetName(ins) != nm->getNdeSetName(ins)) return (1);
  return (0);
}

int exoTestElmBlkNames() {
  int neb = rm->getNumberOfElementBlocks();
  for (int ieb = 1; ieb < neb; ieb++)
    if (rm->getElmBlkName(ieb) != nm->getElmBlkName(ieb)) return (1);
  return (0);
}

// Test implementations
int epTestInputGen(const char *jsonF, const char *ofname, const char *refname) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  for (const auto &prog : inputjson.array_range()) {
    std::unique_ptr<NEM::DRV::NemDriver> nemdrvobj =
        std::unique_ptr<NEM::DRV::NemDriver>(
            NEM::DRV::NemDriver::readJSON(prog));
  }

  return 0;
}

// TEST macros
TEST(EP16, Conversion) {
  EXPECT_EQ(0, exoTestConv(bench1_json, "bctest.g", bench1_ref));
}

TEST(EP16, NodeNumber) { EXPECT_EQ(1, exoTestNodeNum()); }

TEST(EP16, ElementNumber) { EXPECT_EQ(1, exoTestNumElm()); }

TEST(EP16, NodeSetNumber) { EXPECT_EQ(1, exoTestNodeSetNum()); }

TEST(EP16, ElementBlockNumber) { EXPECT_EQ(1, exoTestElmBlkNum()); }

TEST(EP16, TestNodeSetNames) { EXPECT_EQ(0, exoTestNodeSetNames()); }

TEST(EP16, TestElmBlkNames) { EXPECT_EQ(0, exoTestElmBlkNames()); }

TEST(EP16, InputGenerator) {
  EXPECT_EQ(0, epTestInputGen(bench2_json, "bctest.dat", bench2_ref));
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  bench1_json = argv[1];
  bench1_ref = argv[2];
  bench2_json = argv[3];
  bench2_ref = argv[4];

  // run tests
  int res = RUN_ALL_TESTS();

  // clean up
  if (rm) delete rm;
  if (nm) delete nm;

  return res;
}
