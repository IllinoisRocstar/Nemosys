#include <gtest.h>

#include "diffMesh.H"
#include "geoMeshFactory.H"
#include "omegahRefineDriver.H"

// Test cases for NEM::SRV::omegahRefineDriver

const char *refineValueJSON;
const char *refineValueVTU;
const char *refineValueGoldVTU;
const char *refineUniformJSON;
const char *refineUniformVTU;
const char *refineUniformGoldVTU;
const char *refineHexJSON;
const char *refineHexVTU;
const char *refineHexGoldFile;

int genTest(const char *jsonF, const char *newName, const char *goldName) {
  std::string fName(jsonF);
  std::ifstream inputStream(fName);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    return -1;
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  std::unique_ptr<NEM::DRV::omegahRefineDriver> refineDrv(
      NEM::DRV::omegahRefineDriver::readJSON(inputjson));

  auto goldMesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(goldName));
  auto newMesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(newName));

  return NEM::MSH::diffMesh(goldMesh, newMesh);
}

TEST(OmegahRefineDriverTest, RefineUniform) {
  EXPECT_EQ(0,
            genTest(refineUniformJSON, refineUniformVTU, refineUniformGoldVTU));
}

TEST(OmegahRefineDriverTest, RefineValue) {
  EXPECT_EQ(0, genTest(refineValueJSON, refineValueVTU, refineValueGoldVTU));
}

/* TODO: Hex is not yet supported. Requires Omega_h::amr instead of
 *       Omega_h::adapt
TEST(OmegahRefineDriverTest, RefineHexValue) {
  EXPECT_EQ(0, genTest(refineHexJSON, refineHexVTU, refineHexGoldFile));
}
*/

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 10);
  refineValueJSON = argv[1];
  refineValueVTU = argv[2];
  refineValueGoldVTU = argv[3];
  refineUniformJSON = argv[4];
  refineUniformVTU = argv[5];
  refineUniformGoldVTU = argv[6];
  refineHexJSON = argv[7];
  refineHexVTU = argv[8];
  refineHexGoldFile = argv[9];
  return RUN_ALL_TESTS();
}
