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
#ifndef NEMOSYS_COBALT_H_
#define NEMOSYS_COBALT_H_

#include "nemosys_export.h"
#include "Mesh/meshBase.H"

namespace COBALT {

class NEMOSYS_EXPORT cobalt
{
  public:
    cobalt(std::shared_ptr<meshBase> fullMesh,
           const std::string &inFnameVtk,
           const std::string &outFnameCgr,
           const std::string &outFnameCgi);

  public:
    void write() const;
    void writePatchMap(const std::string &mapFile,
                       const std::map<nemId_t, nemId_t> &patchMap) const;
    void writePatchMap(std::ofstream &outputStream,
                       const std::map<nemId_t, nemId_t> &patchMap) const;

  private:
    std::shared_ptr<meshBase> volMeshBase;
    std::shared_ptr<meshBase> surfMeshBase;
    std::string inFnameVtk;
    std::string outFnameCgr;
    std::string outFnameCgi;
};

} // namespace COBALT

#endif  // NEMOSYS_COBALT_H_
