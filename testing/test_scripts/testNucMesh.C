#include <fstream>

#include <gtest.h>
#include <vtkCell.h>
#include <cstdlib>

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
const char *circlesInPolys_test_json;
const char *circlesInPolys_test_REF;
const char *rectangular_array_pattern_test_json;
const char *rectangular_array_pattern_test_REF;
const char *polar_array_pattern_test_json;
const char *polar_array_pattern_test_REF;
const char *hex_array_pattern_test_json;
const char *hex_array_pattern_test_REF;
const char *cartesian_array_test_json;
const char *cartesian_array_test_REF;
const char *mesh_area_conservation_test_json;
const char *mesh_area_conservation_test_REF;
const char *threeD_test_json;
const char *threeD_test_REF;
const char *include_test_json;
const char *include_test_REF;

const double diffGrpTol = 25.0;
const double diffNNdeTol = 2.0;

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

std::string circlesInPolys_test(const char *jsonF) {
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
  std::string ifname = "circlesInPolys_test.msh";
  return ifname;
}

std::string rectangular_array_test(const char *jsonF) {
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
  std::string ifname = "rectangular_array_pattern.msh";
  return ifname;
}

std::string polar_array_test(const char *jsonF) {
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
  std::string ifname = "polar_array_pattern.msh";
  return ifname;
}

std::string hex_array_test(const char *jsonF) {
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
  std::string ifname = "hex_array_pattern_test.msh";
  return ifname;
}

std::string cartesian_array_test(const char *jsonF) {
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
  std::string ifname = "cartesian_array_test.msh";
  return ifname;
}

std::string mesh_area_conservation_test(const char *jsonF) {
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
  std::string ifname = "mesh_area_conservation_test.msh";
  return ifname;
}

std::string threeD_test(const char *jsonF) {
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
  std::string ifname = "3D_test.msh";
  return ifname;
}

std::string include_test(const char *jsonF) {
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
  std::string ifname = "include_test.msh";
  return ifname;
}

// TEST macros
TEST(NucMesh, Simple_Circles_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(2, 0);
  std::vector<int> actualGrp(2, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(simple_circles_test(simple_circles_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(simple_circles_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 2) {
    // actualGrp should equal {120,60}
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Concentric_Circles_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(6, 0);
  std::vector<int> actualGrp(6, 0);

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      concentric_circles_test(concentric_circles_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(concentric_circles_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 6) {
    //  actualGrp should equal {86,90,106,32,20,90};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Concentric_Circles_Test_2) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(6, 0);
  std::vector<int> actualGrp(6, 0);

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      concentric_circles_test_2(concentric_circles_test_2_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(concentric_circles_test_2_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 6) {
    // actualGrp should equal {3498,90,106,32,20,90};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Simple_Polygons_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(9, 0);
  std::vector<int> actualGrp(9, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(simple_polygons_test(simple_polygons_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(simple_polygons_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 9) {
    // actualGrp should equal {18,16,33,12,30,20,33,100,70};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Concentric_Polygons_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(9, 0);
  std::vector<int> actualGrp(9, 0);

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      concentric_polygons_test(concentric_polygons_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(concentric_polygons_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 9) {
    // actualGrp should equal {705,37,45,84,129,62,32,42,292};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, CirclesInPolys_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(4, 0);
  std::vector<int> actualGrp(4, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(circlesInPolys_test(circlesInPolys_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(circlesInPolys_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 4) {
    // actualGrp should equal {512,140,70,140};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Rectangular_Array_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(3, 0);
  std::vector<int> actualGrp(3, 0);

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      rectangular_array_test(rectangular_array_pattern_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(rectangular_array_pattern_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::abs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 3) {
    // actualGrp should equal {3150,800,3156};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Polar_Array_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(5, 0);
  std::vector<int> actualGrp(5, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(polar_array_test(polar_array_pattern_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(polar_array_pattern_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 5) {
    // actualGrp should equal {1187,144,594,2578,160};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Hex_Array_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(3, 0);
  std::vector<int> actualGrp(3, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(hex_array_test(hex_array_pattern_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(hex_array_pattern_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 3) {
    // actualGrp should equal {942,288,360};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Cartesian_Array_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(7, 0);
  std::vector<int> actualGrp(7, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(cartesian_array_test(cartesian_array_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(cartesian_array_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 7) {
    // actualGrp should equal {1604,1848,336,560,280,568,280};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Mesh_Area_Conservation_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(4, 0);
  std::vector<int> actualGrp(4, 0);

  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(
      mesh_area_conservation_test(mesh_area_conservation_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique(mesh_area_conservation_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 4) {
    // actualGrp should equal {12,24,32,9};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, ThreeD_Test) {
  int numNodes1, numNodes2, ret1;
  bool wedge, hex;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(2, 0);
  std::vector<int> actualGrp(2, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(threeD_test(threeD_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_WEDGE) wedge = true;
    if (vtkType == VTK_HEXAHEDRON) hex = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(threeD_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;
  if (diffNNde < diffNNdeTol && wedge == true && hex == true &&
      physGrp_set.size() == 2) {
    // actualGrp should equal {720,360};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

TEST(NucMesh, Include_Test) {
  int numNodes1, numNodes2, ret1;
  bool tris, quads;
  std::set<int> physGrp_set;
  std::vector<int> numElemsInGrp(4, 0);
  std::vector<int> actualGrp(4, 0);

  std::unique_ptr<meshBase> mesh =
      meshBase::CreateUnique(include_test(include_test_json));
  std::vector<double> grpIds(mesh->getNumberOfCells(), 0.0);
  numNodes1 = mesh->getNumberOfPoints();
  mesh->getCellDataArray("PhysGrpId", grpIds);
  for (int iElm = 0; iElm < mesh->getNumberOfCells(); iElm++) {
    VTKCellType vtkType =
        static_cast<VTKCellType>(mesh->getDataSet()->GetCellType(iElm));
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
    // Gather how many elements for each physical group
    numElemsInGrp[grpIds[iElm] - 1] += 1;
  }

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(include_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();
  std::vector<double> grpIdsRef(refMesh->getNumberOfCells(), 0.0);
  refMesh->getCellDataArray("PhysGrpId", grpIdsRef);
  for (int iElm = 0; iElm < refMesh->getNumberOfCells(); iElm++) {
    // Gather how many elements for each physical group
    actualGrp[grpIdsRef[iElm] - 1] += 1;
  }

  double diffNNde = std::fabs(numNodes1 - numNodes2) / numNodes2 * 100.0;

  if (diffNNde < diffNNdeTol && tris == true && quads == true &&
      physGrp_set.size() == 4) {
    // actualGrp should equal {84,64,84,72};
    for (int j = 0; j < numElemsInGrp.size(); ++j) {
      double diffGrp =
          std::fabs(numElemsInGrp[j] - actualGrp[j]) / actualGrp[j] * 100.0;
      if (diffGrp < diffGrpTol)
        ret1 = 0;
      else {
        ret1 = 1;
        std::cout << "Physical Group " << j + 1
                  << " number elements not the same." << std::endl;
        std::cout << "Number of element " << numElemsInGrp[j] << std::endl;
        std::cout << "Number of element actual " << actualGrp[j] << std::endl;
        std::cout << "Difference " << diffGrp << "%" << std::endl;
        break;
      }
    }
  } else
    ret1 = 1;
  ASSERT_EQ(0, ret1);
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 27);
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
  circlesInPolys_test_json = argv[11];
  circlesInPolys_test_REF = argv[12];
  rectangular_array_pattern_test_json = argv[13];
  rectangular_array_pattern_test_REF = argv[14];
  polar_array_pattern_test_json = argv[15];
  polar_array_pattern_test_REF = argv[16];
  hex_array_pattern_test_json = argv[17];
  hex_array_pattern_test_REF = argv[18];
  cartesian_array_test_json = argv[19];
  cartesian_array_test_REF = argv[20];
  mesh_area_conservation_test_json = argv[21];
  mesh_area_conservation_test_REF = argv[22];
  threeD_test_json = argv[23];
  threeD_test_REF = argv[24];
  include_test_json = argv[25];
  include_test_REF = argv[26];

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
  if (!circlesInPolys_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!rectangular_array_pattern_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!polar_array_pattern_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!hex_array_pattern_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!cartesian_array_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!mesh_area_conservation_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!threeD_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }
  if (!include_test_json) {
    std::cerr << "No input file defined" << std::endl;
    throw;
  }

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
