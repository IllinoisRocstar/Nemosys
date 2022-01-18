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
#include <gtest/gtest.h>
#include <Drivers/NemDriver.H>
#include <Mesh/meshBase.H>

const char *bench1_json;
const char *bench1_ref;
const char *bench5_json;
const char *bench5_ref;
const char *bench6_json;
const char *bench6_ref;

int genTest(const char *jsonF, const char *ofname, const char *refname) {
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
    EXPECT_NE(nemdrvobj, nullptr);
    nemdrvobj->execute();
  }

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(refname);
  std::unique_ptr<meshBase> newMesh = meshBase::CreateUnique(ofname);

  return diffMesh(refMesh.get(), newMesh.get());
}

TEST(PNTGen, Bench1) {
  EXPECT_EQ(0, genTest(bench1_json, "bench1_conv.pntmesh", bench1_ref));
}

TEST(PNTGen, Bench5) {
  EXPECT_EQ(0, genTest(bench5_json, "bench5_conv.pntmesh", bench5_ref));
}

TEST(PNTGen, Bench6) {
  EXPECT_EQ(0, genTest(bench6_json, "bench6_conv.pntmesh", bench6_ref));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 7);
  bench1_json = argv[1];
  bench1_ref = argv[2];
  bench5_json = argv[3];
  bench5_ref = argv[4];
  bench6_json = argv[5];
  bench6_ref = argv[6];
  return RUN_ALL_TESTS();
}
