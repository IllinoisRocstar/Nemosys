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
#include <Drivers/NemDriver.H>
#include <gtest/gtest.h>
#include <cctype>
#include <fstream>

const char *vtk2cobaltJSON;
const char *cobalt_out;
const char *cobalt_ref;

TEST(Conversion, ConvertVTKToCobalt) {
  // Test VTK->COBALT conversion driver (does not test the .cgi output)
  // Note quite crude - assumes same ordering of points/elements
  std::ifstream inputStream(vtk2cobaltJSON);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << vtk2cobaltJSON << std::endl;
    exit(1);
  }
  jsoncons::json inputjson;
  inputStream >> inputjson;
  NEM::DRV::NemDriver::readJSON(inputjson)->execute();

  std::ifstream resultFile(cobalt_out);
  EXPECT_TRUE(resultFile.good());
  std::ifstream reference(cobalt_ref);
  EXPECT_TRUE(reference.good());
  {
    std::string lineRef, lineResult;
    while (std::getline(reference, lineRef)) {
      EXPECT_TRUE(std::getline(resultFile, lineResult));
      std::stringstream refStream(lineRef);
      std::stringstream resultStream(lineResult);
      if (lineRef.empty()) {
        EXPECT_TRUE(lineResult.empty());
        EXPECT_TRUE(reference.eof());
        EXPECT_TRUE(resultFile.eof());
      } else if (std::isspace(lineRef.at(0))) {
        // Line beginning with whitespace => coordinates
        EXPECT_TRUE(std::isspace(lineResult.at(0)));
        while (true) {
          float refVal;
          refStream >> refVal;
          float resultVal;
          resultStream >> resultVal;
          if (!refStream) {
            EXPECT_FALSE(resultStream);
            break;
          } else {
            EXPECT_TRUE(resultStream);
          }
          EXPECT_FLOAT_EQ(refVal, resultVal);
        }
      } else {
        // No white space => indices/counts
        EXPECT_FALSE(std::isspace(lineResult.at(0)));
        while (true) {
          long long int refVal;
          refStream >> refVal;
          long long int resultVal;
          resultStream >> resultVal;
          if (!refStream) {
            EXPECT_FALSE(resultStream);
            break;
          } else {
            EXPECT_TRUE(resultStream);
          }
          EXPECT_EQ(refVal, resultVal);
        }
      }
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 4);
  vtk2cobaltJSON = argv[1];
  cobalt_out = argv[2];
  cobalt_ref = argv[3];
  return RUN_ALL_TESTS();
}
