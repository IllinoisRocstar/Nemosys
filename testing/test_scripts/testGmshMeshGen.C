#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include <fstream>

#include <gtest.h>
#include <vtkCell.h>

#include "Drivers/MeshGen/GmshMeshGenDriver.H"
#include "Drivers/NemDriver.H"
#include "geoMeshFactory.H"
#include "meshBase.H"

const char *box_test_json;
const char *pitz_daily_test_json;
const char *pitz_daily_test_REF;
const char *box_test_REF;

// Test implementations
std::string box_test(const char *jsonF) {
  std::cout << "Starting box_test method\n" << std::endl;
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }
  jsoncons::json inputjson;
  inputStream >> inputjson;

  for (const auto &prog : inputjson.array_range()) {
    std::cout << "Reading JSON array" << std::endl;
    auto nemdrvobj = NEM::DRV::NemDriver::readJSON(prog);
    EXPECT_NE(dynamic_cast<NEM::DRV::GmshMeshGenDriver *>(nemdrvobj.get()), nullptr);
    nemdrvobj->execute();
  }

  std::string ifname = "./box_test.vtu";
  std::cout << "Output mesh file from box_test " << ifname << std::endl;
  return ifname;
}

TEST(gmshMeshGenTest, pitz_daily_Test) {
  std::cout << "Running flow box test" << std::endl;
  std::string fname(pitz_daily_test_json);
  std::ifstream inputStream(fname);
  jsoncons::json inputjson;
  inputStream >> inputjson;
  for (const auto &prog : inputjson.array_range()) {
    auto nemdrvobj = NEM::DRV::NemDriver::readJSON(prog);
    EXPECT_NE(nemdrvobj, nullptr);
    nemdrvobj->execute();
  }

  auto genPDMesh = std::shared_ptr<NEM::MSH::geoMeshBase>(NEM::MSH::Read("pitzdaily.msh"));
  auto refPDMesh = std::shared_ptr<NEM::MSH::geoMeshBase>(NEM::MSH::Read(pitz_daily_test_REF));
  EXPECT_NEAR(refPDMesh->getNumberOfPoints(),
              genPDMesh->getNumberOfPoints(),
              .1*refPDMesh->getNumberOfPoints());
  EXPECT_NEAR(refPDMesh->getNumberOfCells(),
              genPDMesh->getNumberOfCells(),
              .1*refPDMesh->getNumberOfCells());
}
// TEST macros
TEST(gmshMeshGenTest, Box_Test) {

  int newNodes = 0;
  int refNodes = 0;
  int ret1 = 0;

  std::cout << "Running box_test test" << std::endl;
  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(box_test(box_test_json));
  newNodes = mesh->getNumberOfPoints();
  std::cout << "newNodes = " << newNodes << std::endl;

  std::cout << box_test_REF << std::endl;
  std::cout << "Reading reference mesh file" << std::endl;
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(box_test_REF);
  refNodes = refMesh->getNumberOfPoints();
  std::cout << "refNodes = " << refNodes << std::endl;

#ifdef _WIN32
  int d = 50;
#else
  int d = 100;
#endif

  if (refNodes - refNodes / d <= newNodes &&
      newNodes <= refNodes + refNodes / d)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  box_test_json = argv[1];
  box_test_REF = argv[2];

  pitz_daily_test_json = argv[3];
  pitz_daily_test_REF = argv[4];

  if (!box_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
