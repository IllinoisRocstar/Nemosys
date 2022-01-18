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
#ifndef NEMOSYS_MESHSTITCHER_H_
#define NEMOSYS_MESHSTITCHER_H_

//std
#include <vector>
#include <string>
#include <memory>

#include "nemosys_export.h"

class cgnsAnalyzer;
class rocstarCgns;
class meshBase;

/* This class stitches together cgns grids and solution data. It exports the 
   stitched mesh with data to VTK data structures. The 1st component of 
   the member partitions, partitions[0], contains the stitched cgns object
   after instantiation of meshStitcher (when stitching vol/fluid files). 
    cgObj is the stitched mesh when stitching surf/ifluid files */
class NEMOSYS_EXPORT meshStitcher {
  public:
    // stitch fluid/ifluid files depending on surf flag
    // surface stitching uses rocstarCgns, while vol uses cgnsAnalyzer 
    meshStitcher(std::vector<std::string> cgFileNames, bool surf);

    ~meshStitcher() = default;

    // disabling copy and copy-assignment operators
    meshStitcher(const meshStitcher &that) = delete;
    meshStitcher &operator=(const meshStitcher &that) = delete;

    std::shared_ptr<cgnsAnalyzer> getStitchedCGNS() const;
    std::shared_ptr<meshBase> getStitchedMB() const;

  private:
    // names of cgns files
    const std::vector<std::string> cgFileNames;
    // vector of loaded cgns files
    std::vector<std::shared_ptr<cgnsAnalyzer>> partitions;
    // stitched mesh in Nemosys format  
    std::shared_ptr<meshBase> stitchedMesh;
    // stitched mesh in cgns fomrat
    std::shared_ptr<rocstarCgns> cgObj;

  private:
    // initialize series of cgns objects
    void initSurfCgObj();
    void initVolCgObj();
};

#endif  // NEMOSYS_MESHSTITCHER_H_
