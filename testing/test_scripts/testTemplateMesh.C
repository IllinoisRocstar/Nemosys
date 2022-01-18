/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#include <fstream>
#include <string>
#include <gtest/gtest.h>
#include <vtkCell.h>

#include <Drivers/NemDriver.H>
#include <Mesh/meshBase.H>

const char *test_spiral_tape_json;
const char *spiral_test_REF;

// Test implementations
std::string spiral_tape_test(const std::string &jsonF) {
  std::ifstream inputStream(jsonF);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }
  jsoncons::json inputjson;
  inputStream >> inputjson;

  auto nemdrvobj = NEM::DRV::NemDriver::readJSON(inputjson);
  nemdrvobj->execute();

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
    if (vtkType == VTK_TRIANGLE) tris = true;
    if (vtkType == VTK_QUAD) quads = true;
    if (vtkType == VTK_HEXAHEDRON) hexs = true;
    physGrp_set.insert(static_cast<int>(grpIds[iElm]));
  }

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(spiral_test_REF);
  numNodes2 = refMesh->getNumberOfPoints();

  if (std::abs(numNodes1 - numNodes2) < 0.05 * numNodes2 && tris == false &&
      quads == true && hexs == true && physGrp_set.size() == 5)
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
