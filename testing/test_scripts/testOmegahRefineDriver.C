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
#include <Drivers/Refine/OmegahRefineDriver.H>
#include <Mesh/diffMesh.H>
#include <Mesh/geoMeshFactory.H>

// Test cases for NEM::SRV::omegahRefineDriver

const char *refineValueJSON;
const char *refineValueVTU;
const char *refineValueGoldVTU;
const char *refineUniformJSON;
const char *refineUniformVTU;
const char *refineUniformGoldVTU;
const char *refineHexJSON;
const char *refineHexVTU;
const char *refineHexGoldFile;

int genTest(const char *jsonF, const char *newName, const char *goldName) {
  std::string fName(jsonF);
  std::ifstream inputStream(fName);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    return -1;
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  auto refineDrv = NEM::DRV::NemDriver::readJSON(inputjson);
  EXPECT_NE(dynamic_cast<NEM::DRV::OmegahRefineDriver *>(refineDrv.get()),
            nullptr);
  refineDrv->execute();

  auto goldMesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(goldName));
  auto newMesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(newName));

  return NEM::MSH::diffMesh(goldMesh, newMesh);
}

TEST(OmegahRefineDriverTest, RefineUniform) {
  EXPECT_EQ(0,
            genTest(refineUniformJSON, refineUniformVTU, refineUniformGoldVTU));
}

TEST(OmegahRefineDriverTest, RefineValue) {
  EXPECT_EQ(0, genTest(refineValueJSON, refineValueVTU, refineValueGoldVTU));
}

/* TODO: Hex is not yet supported. Requires Omega_h::amr instead of
 *       Omega_h::adapt
TEST(OmegahRefineDriverTest, RefineHexValue) {
  EXPECT_EQ(0, genTest(refineHexJSON, refineHexVTU, refineHexGoldFile));
}
*/

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 10);
  refineValueJSON = argv[1];
  refineValueVTU = argv[2];
  refineValueGoldVTU = argv[3];
  refineUniformJSON = argv[4];
  refineUniformVTU = argv[5];
  refineUniformGoldVTU = argv[6];
  refineHexJSON = argv[7];
  refineHexVTU = argv[8];
  refineHexGoldFile = argv[9];
  return RUN_ALL_TESTS();
}
