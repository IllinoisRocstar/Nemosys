#include <fstream>

#include <gtest.h>
#include <vtkCell.h>

#include "NemDriver.H"
#include "meshBase.H"

const char *simple_circles_test_json;
const char *simple_circles_test_REF;
const char *concentric_circles_test_json;
const char *concentric_circles_test_REF;
const char *concentric_circles_test_2_json;
const char *concentric_circles_test_2_REF;
const char *simple_polygons_test_json;
const char *simple_polygons_test_REF;
const char *concentric_polygons_test_json;
const char *concentric_polygons_test_REF;

// Test implementations
std::string simple_circles_test(const char *jsonF) {
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
  std::string ifname = "simple_circles_test.msh";
  return ifname;
}

std::string concentric_circles_test(const char *jsonF) {
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
  std::string ifname = "concentric_circles_test.msh";
  return ifname;
}

std::string concentric_circles_test_2(const char *jsonF) {
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
  std::string ifname = "concentric_circles_test_2.msh";
  return ifname;
}

std::string simple_polygons_test(const char *jsonF) {
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
  std::string ifname = "simple_polygons_test.msh";
  return ifname;
}

std::string concentric_polygons_test(const char *jsonF) {
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
  std::string ifname = "concentric_polygons_test.msh";
  return ifname;
}

// TEST macros
TEST(NucMesh, Simple_Circles_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(simple_circles_test(simple_circles_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE)
      tris = true;
    if (vtkType == VTK_QUAD)
      quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(simple_circles_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();

  if (numNodes1 == numNodes2 && tris == true && quads == true &&
      physGrp_set.size() == 2)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Concentric_Circles_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      concentric_circles_test(concentric_circles_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE)
      tris = true;
    if (vtkType == VTK_QUAD)
      quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(concentric_circles_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();

  if (numNodes1 == numNodes2 && tris == true && quads == true &&
      physGrp_set.size() == 6)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Concentric_Circles_Test_2) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      concentric_circles_test_2(concentric_circles_test_2_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE)
      tris = true;
    if (vtkType == VTK_QUAD)
      quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(concentric_circles_test_2_REF);
  numNodes2 = refMesh->getNumberOfPoints();

  if (numNodes1 == numNodes2 && tris == true && quads == true &&
      physGrp_set.size() == 6)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Simple_Polygons_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(simple_polygons_test(simple_polygons_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE)
      tris = true;
    if (vtkType == VTK_QUAD)
      quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(simple_polygons_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();

  if (numNodes1 == numNodes2 && tris == true && quads == true &&
      physGrp_set.size() == 9)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Concentric_Polygons_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      concentric_polygons_test(concentric_polygons_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE)
      tris = true;
    if (vtkType == VTK_QUAD)
      quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(concentric_polygons_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();

  if (numNodes1 == numNodes2 && tris == true && quads == true &&
      physGrp_set.size() == 9)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 11);
  simple_circles_test_json = argv[1];
  simple_circles_test_REF = argv[2];
  concentric_circles_test_json = argv[3];
  concentric_circles_test_REF = argv[4];
  concentric_circles_test_2_json = argv[5];
  concentric_circles_test_2_REF = argv[6];
  simple_polygons_test_json = argv[7];
  simple_polygons_test_REF = argv[8];
  concentric_polygons_test_json = argv[9];
  concentric_polygons_test_REF = argv[10];

  if (!simple_circles_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!concentric_circles_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!concentric_circles_test_2_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!simple_polygons_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!concentric_polygons_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
