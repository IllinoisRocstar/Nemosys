#include <fstream>
#include <stdio.h>

#include <gtest.h>
#include <vtkCell.h>

#include "NemDriver.H"
#include "meshBase.H"

const char *test_spiral_tape_json;
const char *spiral_test_REF;

// Test implementations
std::string spiral_tape_test(const char *jsonF) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }
  jsoncons::json inputjson;
  inputStream >> inputjson;

  if (inputjson.contains("Template Name")) {
    std::unique_ptr<NemDriver> nemdrvobj =
        std::unique_ptr<NemDriver>(NemDriver::readJSON(inputjson));
  }
  std::string ifname = "spiral.msh";
  return ifname;
}


// TEST macros
TEST(TemplateMesh, Spiral_Tape_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads, hexs;
  std::set<int> physGrp_set;

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(spiral_tape_test(test_spiral_tape_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  tris = false;
  quads = false;
  hexs = false;
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE)
      tris = true;
    if (vtkType == VTK_QUAD)
      quads = true;
    if (vtkType == VTK_HEXAHEDRON)
      hexs = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(spiral_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();

  if (numNodes1 == numNodes2 && tris == false && quads == true
      && hexs == true && physGrp_set.size() == 5)
    ret1 = 0;
  else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 3);
  test_spiral_tape_json = argv[1];
  spiral_test_REF = argv[2];

  if (!test_spiral_tape_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
