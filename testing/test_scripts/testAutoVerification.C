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
#include "omp.h"

#include <Drivers/AutoVerificationDriver.H>
#include <SolutionVerification/OrderOfAccuracy.H>

const char *coarse;
const char *fine;
const char *finer;
const char *rich_ref;
const char *jsonfile;

// TODO Test Description

// TODO set EQUAL to NEAR
// TODO set num threads to 2 for test

class OrderOfAccuracyTest : public ::testing::Test {
 protected:
  OrderOfAccuracyTest() {
    f3 = meshBase::CreateUnique(coarse);
    f2 = meshBase::CreateUnique(fine);
    f1 = meshBase::CreateUnique(finer);
    oac = new OrderOfAccuracy(f3.get(), f2.get(), f1.get(), arrayIDs);
  }

  virtual ~OrderOfAccuracyTest() {}

  virtual void SetUp() {}

  virtual void TearDown() {}

  // TODO get rid of unique pointers -> convert to shared
  std::unique_ptr<meshBase> f1;
  std::unique_ptr<meshBase> f2;
  std::unique_ptr<meshBase> f3;
  const std::vector<int> arrayIDs = {0};
  const std::vector<std::vector<double>> res_ref = {
      {0.32581996578105021, 0.19338485350409632, 0.32030546671960047}};
  const std::vector<std::vector<double>> oac_ref = {
      {2.1119174483287679, 1.9704754800492528, 2.0774696352142827}};
  const std::vector<std::vector<double>> asymp_ref = {
      {1.0200502831373015, 1.0390518032820475, 1.031263791985966}};
  OrderOfAccuracy *oac;
};

// tests Richardson extrapolation
TEST_F(OrderOfAccuracyTest, ComputeRichardsonExtrap) {
  oac->computeRichardsonExtrapolation();
  std::unique_ptr<meshBase> f1_ref = meshBase::CreateUnique(rich_ref);
  EXPECT_EQ(0, diffMesh(f1.get(), f1_ref.get()));
}

// tests if correct resolution desired error is returned
TEST_F(OrderOfAccuracyTest, ComputeResolution) {
  std::vector<std::vector<double>> res = oac->computeResolution(.005);
  for (int j = 0; j < res[0].size(); ++j) {
    EXPECT_NEAR(res_ref[0][j], res[0][j], 1e-12);
  }
}

// tests if correct oac is returned
TEST_F(OrderOfAccuracyTest, ComputeOrderOfAccuracy) {
  std::vector<std::vector<double>> order = oac->computeOrderOfAccuracy();
  for (int j = 0; j < order[0].size(); ++j) {
    EXPECT_NEAR(order[0][j], oac_ref[0][j], 1e-12);
  }
}

// tests if correct asymptotic ratios are returned
TEST_F(OrderOfAccuracyTest, CheckAsymptoticRange) {
  std::vector<std::vector<double>> asymp = oac->checkAsymptoticRange();
  for (int j = 0; j < asymp[0].size(); ++j) {
    EXPECT_NEAR(asymp[0][j], asymp_ref[0][j], 1e-12);
  }
}

TEST(AutoVerificationDriverTest, DriverTest) {
  std::ifstream inputStream(jsonfile);
  jsoncons::json inputjson;
  inputStream >> inputjson;
  auto driver = NEM::DRV::NemDriver::readJSON(inputjson);
  EXPECT_NE(driver, nullptr);
  driver->execute();
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 6);
  finer = argv[1];
  fine = argv[2];
  coarse = argv[3];
  rich_ref = argv[4];
  jsonfile = argv[5];
  cout << "Max number of threads available " << omp_get_max_threads() << std::endl;
  return RUN_ALL_TESTS();
}
