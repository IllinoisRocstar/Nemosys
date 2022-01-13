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
#include <gtest/gtest.h>
#include <stdlib.h>
#include <vtkPointData.h>
#include <Drivers/ProteusDriver.H>
#include <Mesh/meshBase.H>
#include <Mesh/exoMesh.H>

const char *arg_jsonFile1;
std::string arg_vtuFile1;
std::string arg_vtuFileGold1;
std::string arg_exoFile1;
std::string arg_exoFileGold1;

int convert(const char *jsonF) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  for (const auto &prog : inputjson.array_range()) {
    auto nemdrvobj = NEM::DRV::NemDriver::readJSON(prog);
    EXPECT_NE(dynamic_cast<NEM::DRV::ProteusDriver *>(nemdrvobj.get()),
              nullptr);
    nemdrvobj->execute();
  }

  return 0;
}

int test_diffExoMesh(const std::string &fName1, const std::string &fName2) {
  NEM::MSH::EXOMesh::exoMesh em1;
  NEM::MSH::EXOMesh::exoMesh em2;

  em1.read(fName1);
  em2.read(fName2);

  if (em1.getNumberOfNodes() != em2.getNumberOfNodes()) {
    std::cout << "Meshes don't have the same number of nodes" << std::endl;
    return 1;
  }
  if (em1.getNumberOfElements() != em2.getNumberOfElements()) {
    std::cout << "Meshes don't have the same number of elements" << std::endl;
    return 1;
  }
  if (em1.getNumberOfNodeSets() != em2.getNumberOfNodeSets()) {
    std::cout << "Meshes don't have the same number of node sets" << std::endl;
    return 1;
  }
  if (em1.getNumberOfElementBlocks() != em2.getNumberOfElementBlocks()) {
    std::cout << "Meshes don't have the same number of element blocks"
              << std::endl;
    return 1;
  }
  if (em1.getNumberOfSideSets() != em2.getNumberOfSideSets()) {
    std::cout << "Meshes don't have the same number of side sets" << std::endl;
    return 1;
  }
  std::cout << "Meshes are the same" << std::endl;
  return 0;
}

int test_diffVtkMesh(const std::string &fName1, const std::string &fName2,
                     double tol) {
  meshBase *mesh1 = meshBase::Create(fName1);
  meshBase *mesh2 = meshBase::Create(fName2);

  if (mesh1->getNumberOfPoints() != mesh2->getNumberOfPoints() ||
      mesh1->getNumberOfCells() != mesh2->getNumberOfCells()) {
    std::cout << "Meshes don't have the same number of points or cells"
              << std::endl;
    return 1;
  }

  for (int i = 0; i < mesh1->getNumberOfPoints(); ++i) {
    std::vector<double> coord1 = mesh1->getPoint(i);
    std::vector<double> coord2 = mesh2->getPoint(i);
    for (int j = 0; j < 3; ++j) {
      if (std::fabs(coord1[j] - coord2[j]) > tol) {
        std::cout << "Meshes differ in point coordinates" << std::endl;
        return 1;
      }
    }
  }

  for (int i = 0; i < mesh1->getNumberOfCells(); ++i) {
    std::vector<std::vector<double>> cell1 = mesh1->getCellVec(i);
    std::vector<std::vector<double>> cell2 = mesh2->getCellVec(i);
    if (cell1.size() != cell2.size()) {
      std::cout << "Meshes differ in cells" << std::endl;
      return 1;
      ;
    }
    for (int j = 0; j < cell1.size(); ++j) {
      for (int k = 0; k < 3; ++k) {
        if (std::fabs(cell1[j][k] - cell2[j][k]) > tol) {
          std::cout << "Meshes differ in cells" << std::endl;
          return 1;
        }
      }
    }
  }

  vtkSmartPointer<vtkPointData> pd1 = vtkSmartPointer<vtkPointData>::New();
  pd1 = mesh1->getDataSet()->GetPointData();
  vtkSmartPointer<vtkPointData> pd2 = vtkSmartPointer<vtkPointData>::New();
  pd2 = mesh2->getDataSet()->GetPointData();
  int numArr1 = pd1->GetNumberOfArrays();
  int numArr2 = pd2->GetNumberOfArrays();

  if (numArr1 != numArr2) {
    std::cout << "Meshes have different numbers of point data" << std::endl;
    return 1;
  }

  for (int i = 0; i < numArr1; ++i) {
    vtkDataArray *da1 = pd1->GetArray(i);
    vtkDataArray *da2 = pd2->GetArray(i);
    int numComponent = da1->GetNumberOfComponents();
    for (int j = 0; j < mesh1->getNumberOfPoints(); ++j) {
      std::vector<double> comps1(numComponent);
      std::vector<double> comps2(numComponent);
      da1->GetTuple(j, comps1.data());
      da2->GetTuple(j, comps2.data());
      for (int k = 0; k < numComponent; ++k) {
        if (std::fabs(comps1[k] - comps2[k]) > tol) {
          std::cout << "Meshes differ in point data values" << std::endl;
          return 1;
        }
      }
    }
  }
  if (mesh1)
    delete mesh1;
  if (mesh2)
    delete mesh2;
  std::cout << "Meshes are the same" << std::endl;
  return 0;
}

TEST(ProteusDriverTest, convert1) { EXPECT_EQ(0, convert(arg_jsonFile1)); }
TEST(exoMesh, diffExo1) {
  EXPECT_EQ(0, test_diffExoMesh(arg_exoFile1, arg_exoFileGold1));
}
TEST(vtkMesh, diffVtk1) {
  EXPECT_EQ(0, test_diffVtkMesh(arg_vtuFile1, arg_vtuFileGold1, 1e-8));
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 6);

  arg_jsonFile1 = argv[1];
  arg_vtuFile1 = argv[2];
  arg_vtuFileGold1 = argv[3];
  arg_exoFile1 = argv[4];
  arg_exoFileGold1 = argv[5];

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
