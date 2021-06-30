#include <fstream>

#include <gtest.h>

#include "Drivers/NemDriver.H"
#include "meshBase.H"

const char *defaultjson;
const char *defaultRef;
const char *unifjson;
const char *unifRef;
const char *geomjson;
const char *geomRef;

int genTest(const char *jsonF, const char *ofname, const char *refname) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << defaultjson << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  for (const auto &prog : inputjson.array_range()) {
    auto nemdrvobj = NEM::DRV::NemDriver::readJSON(prog);
    EXPECT_NE(nemdrvobj, nullptr);
    nemdrvobj->execute();
  }

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(refname);
  std::unique_ptr<meshBase> newMesh = meshBase::CreateUnique(ofname);

  // return diffMesh(refMesh.get(),newMesh.get());

  // Due to Netgen algorithm changing slightly between different versions,
  // we are not going to compare the meshes point-by-point, cell-by-cell.
  // The test will check the number of cells and points and ensure they
  // are within a 0.5% tolerance.

#ifdef _WIN32
  nemId_t divisor = 39;  // 2.56 % = 1 / 39
#else
  nemId_t divisor = 200;  // 0.5% = 1 / 200
#endif

  nemId_t refpoints = refMesh->getNumberOfPoints();
  nemId_t refcells = refMesh->getNumberOfCells();
  nemId_t newpoints = newMesh->getNumberOfPoints();
  nemId_t newcells = newMesh->getNumberOfCells();

  std::cout << "refpoints = " << refpoints << std::endl;
  std::cout << "refcells  = " << refcells << std::endl;
  std::cout << "newPoints = " << newpoints << std::endl;
  std::cout << "newcells  = " << newcells << std::endl;
  std::cout << "Percent difference in points = "
            << 100.0 * (std::abs(static_cast<double>(refpoints) -
                                 static_cast<double>(newpoints)) /
                        static_cast<double>(refpoints))
            << std::endl;
  std::cout << "Percent difference in cells = "
            << 100.0 * (std::abs(static_cast<double>(refcells) -
                                 static_cast<double>(newcells)) /
                        static_cast<double>(refcells))
            << std::endl;

  return !(refpoints - refpoints / divisor <= newpoints &&
           newpoints <= refpoints + refpoints / divisor &&
           refcells - refcells / divisor <= newcells &&
           newcells <= refcells + refcells / divisor);
}

TEST(MeshGen, defaultGen) {
  EXPECT_EQ(0, genTest(defaultjson, "hinge.vtu", defaultRef));
}

TEST(MeshGen, uniformRefinementGen) {
  EXPECT_EQ(0, genTest(unifjson, "hingeUnif.vtu", unifRef));
}

TEST(MeshGen, geometryRefinementGen) {
  EXPECT_EQ(0, genTest(geomjson, "hingeGeom.vtu", geomRef));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 7);
  defaultjson = argv[1];
  defaultRef = argv[2];
  unifjson = argv[3];
  unifRef = argv[4];
  geomjson = argv[5];
  geomRef = argv[6];
  return RUN_ALL_TESTS();
}
