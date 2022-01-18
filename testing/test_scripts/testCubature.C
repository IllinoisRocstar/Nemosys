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
#include <Integration/Cubature.H>
#include <Mesh/meshBase.H>

const char *nodeMesh;
const char *refGauss;
const char *refGauss1;
const char *refInt;
const char *hexmesh;

// The fixture for testing class orthoPoly. From google test primer.
class CubatureTest : public ::testing::Test {
 protected:
  // You can remove any or all of the following functions if its body
  // is empty.

  CubatureTest() {
    // mesh = meshBase::Create(nodeMesh);
    // cuby = new GaussCubature(mesh,arrayIDs);

    mesh = meshBase::CreateShared(nodeMesh);
    cuby = GaussCubature::CreateShared(mesh->getDataSet(), arrayIDs);
  }

  virtual ~CubatureTest() {}

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:
  virtual void SetUp() {}

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects declared here can be used by all tests in the test case for
  // orthoPoly.
  std::shared_ptr<meshBase> mesh;
  std::shared_ptr<GaussCubature> cuby;
  const std::vector<int> arrayIDs = {0, 2, 3, 7};
};

TEST(CubatureContructors, ConstructWithoutArray) {
  // load reference node mesh
  std::unique_ptr<meshBase> mesh1 = meshBase::CreateUnique(nodeMesh);
  // generate gauss point mesh
  std::unique_ptr<GaussCubature> cubeObj =
      GaussCubature::CreateUnique(mesh1->getDataSet());
  cubeObj->writeGaussMesh("gaussTestNoData.vtp");
  // load generated mesh
  std::unique_ptr<meshBase> gaussMesh =
      meshBase::CreateUnique("gaussTestNoData.vtp");
  // load reference gauss mesh
  std::unique_ptr<meshBase> refGaussMesh = meshBase::CreateUnique(refGauss1);
  EXPECT_EQ(0, diffMesh(gaussMesh.get(), refGaussMesh.get()));
}

// tests construction with array, writeout, loading and interpolation
TEST_F(CubatureTest, InterpolateToGauss) {
  cuby->writeGaussMesh("gaussTest.vtp");
  std::unique_ptr<meshBase> gaussMesh = meshBase::CreateUnique("gaussTest.vtp");
  std::unique_ptr<meshBase> refGaussMesh = meshBase::CreateUnique(refGauss);
  EXPECT_EQ(0, diffMesh(gaussMesh.get(), refGaussMesh.get()));
}

// test integration methods
TEST_F(CubatureTest, IntegrateOverAllCells) {
  std::vector<std::vector<double>> totalIntegralData(
      cuby->integrateOverAllCells());
  std::vector<std::vector<double>> refTotalIntegralData = {
      {28385.23446875294},
      {12718.489406555998},
      {10811985.874199742},
      {0.00042734000024760349, -5.7540334078832547e-08,
       -5.9937717104360543e-08}};

  for (int i = 0; i < totalIntegralData.size(); ++i) {
    for (int j = 0; j < totalIntegralData[i].size(); ++j) {
      EXPECT_NEAR(totalIntegralData[i][j], refTotalIntegralData[i][j],1.E-5);
    }
  }

  std::unique_ptr<meshBase> integralMesh = meshBase::CreateUnique(refInt);
  std::unique_ptr<meshBase> cubatureMesh =
      meshBase::CreateUnique(cuby->getDataSet(), "Cubature");
  EXPECT_EQ(0, diffMesh(cubatureMesh.get(), integralMesh.get()));
}

TEST(CubatureHexTest, integrateOverHex) {
  std::unique_ptr<meshBase> hex = meshBase::CreateUnique(hexmesh);
  std::vector<int> arrayIDs = {0};
  std::unique_ptr<GaussCubature> cubeObj =
      GaussCubature::CreateUnique(hex->getDataSet(), arrayIDs);
  std::vector<std::vector<double>> totalIntegralData(
      cubeObj->integrateOverAllCells());
  EXPECT_DOUBLE_EQ(totalIntegralData[0][0], 0.01780325297629036);
/*
  std::cout << setprecision(15) << totalIntegralData[0][0] << std::endl;
  double p = 0.577350269189626;
  std::vector<double> gp1 = {-p, -p, -p};
  std::vector<double> gp2 = {p, -p, -p};
  std::vector<double> gp3 = {-p, p, -p};
  std::vector<double> gp4 = {p, p, -p};
  std::vector<double> gp5 = {-p, -p, p};
  std::vector<double> gp6 = {p, -p, p};
  std::vector<double> gp7 = {-p, p, p};
  std::vector<double> gp8 = {p, p, p};

  std::vector<std::vector<double>> gp(8);
  gp[0] = gp1;
  gp[1] = gp2;
  gp[2] = gp3;
  gp[3] = gp4;
  gp[4] = gp5;
  gp[5] = gp6;
  gp[6] = gp7;
  gp[7] = gp8;

  vtkSmartPointer<vtkGenericCell> genCell =
      vtkSmartPointer<vtkGenericCell>::New();
  hex->getDataSet()->GetCell(0, genCell);
  double pcoords[3];
  double weights[8];
  double minDist2;
  int subId;
  for (int i = 0; i < 8; ++i) {
    genCell->EvaluatePosition(gp[i].data(), NULL, subId, pcoords, minDist2,
                              weights);
    for (int j = 0; j < 8; ++j) {
      std::cout << std::setprecision(15) << weights[j] << ", ";
    }
    std::cout << std::endl;
  }
*/
}

/*
double integrand(const std::vector<double> &coord) {
  return sin(coord[0]) + coord[1] * coord[1] + cos(coord[2]);
}

TEST(IntegrationTest, writeNewCube) {
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique("single-hex.vtk");
  vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
  da->SetName("integrand");
  da->SetNumberOfComponents(1);
  da->SetNumberOfTuples(mesh->getNumberOfPoints());

  for (int i = 0; i < mesh->getNumberOfPoints(); ++i) {
    double val[1];
    val[0] = integrand(mesh->getPoint(i));
    da->InsertTuple(i, val);
  }
  mesh->getDataSet()->GetPointData()->AddArray(da);
  mesh->write();
}
*/

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 6);
  nodeMesh = argv[1];
  refGauss = argv[2];
  refGauss1 = argv[3];
  refInt = argv[4];
  hexmesh = argv[5];
  return RUN_ALL_TESTS();
}

/*
double integrand(const std::vector<double> &coord) {
  return sin(coord[0]) + coord[1] * coord[1] + cos(coord[2]);
}

double integrand1(const std::vector<double> &coord) {
  return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
}
*/
