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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: " << argv[0] << " input.json" << std::endl;
    exit(1);
  }

  std::string fname = argv[1];
  auto ext_begin = fname.find_last_of('.');
  if (ext_begin == std::string::npos ||
      fname.compare(ext_begin, fname.length() - ext_begin, ".json") != 0) {
    std::cerr << "Warning: " << fname << " does not end in .json\n";
  }

  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << fname << std::endl;
    std::exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  if (inputjson.is_array())
    for (const auto &prog : inputjson.array_range()) {
      auto nemdrvobj = NEM::DRV::NemDriver::readJSON(prog);
      nemdrvobj->execute();
    }
  else {
    auto nemdrvobj = NEM::DRV::NemDriver::readJSON(inputjson);
    nemdrvobj->execute();
  }

  return 0;
}

// TODO: also overload each driver constructor with one that takes pointer to
//       resulting mesh with that, we can maybe do more things in memory.
//       Add Netgen uniform option to Refinement program (i.e. Refinement Method
//       = "uniformMAd" || "uniformNG")
