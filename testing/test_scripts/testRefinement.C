#include <RefineDriver.H>
#include <gtest.h>

const char *refineValueJSON;
const char *refineValueVTU;
const char *refineValueGoldVTU;
const char *refineUniformJSON;
const char *refineUniformVTU;
const char *refineUniformGoldVTU;

int genTest(const char *jsonF, const char *newName, const char *goldName) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  std::unique_ptr<RefineDriver> refineDriver =
      std::unique_ptr<RefineDriver>(RefineDriver::readJSON(inputjson));

  std::unique_ptr<meshBase> goldMesh = meshBase::CreateUnique(goldName);
  std::unique_ptr<meshBase> newMesh = meshBase::CreateUnique(newName);

  // return 0;
  return diffMesh(goldMesh.get(), newMesh.get());
}

TEST(RefinementDriverTest, RefineValue) {
  EXPECT_EQ(0, genTest(refineValueJSON, refineValueVTU, refineValueGoldVTU));
}

TEST(RefinementDriverTest, RefineUniform) {
  EXPECT_EQ(0,
            genTest(refineUniformJSON, refineUniformVTU, refineUniformGoldVTU));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 7);
  refineValueJSON = argv[1];
  refineValueVTU = argv[2];
  refineValueGoldVTU = argv[3];
  refineUniformJSON = argv[4];
  refineUniformVTU = argv[5];
  refineUniformGoldVTU = argv[6];
  return RUN_ALL_TESTS();
}
