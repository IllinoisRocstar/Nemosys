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

#include <InputGeneration/ep16Post.H>
#include <vtkSTLReader.h>
#include <string>

// data loading
class TestDriver : public testing::Test {
 protected:
  jsoncons::json inp;

  virtual void SetUp() {
    std::ifstream inputStream;
    inputStream.open("ep16post.json");
    inputStream >> inp;
  }
};

TEST_F(TestDriver, EndToEnd) {
  int ret = 0;
  NEM::EPC::ep16Post *ep =
      NEM::EPC::ep16Post::readJSON(inp.at("Input Generation Options"), ret);
  ASSERT_EQ(0, ret);
  delete ep;
  vtkNew<vtkSTLReader> stlReader;
  for (int i = 0; i < 4; ++i) {
    std::string stlFileName{"clust_" + std::to_string(i) + ".stl"};
    stlReader->SetFileName(stlFileName.c_str());
    stlReader->Update();
    ASSERT_NE(stlReader->GetOutput(), nullptr);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
