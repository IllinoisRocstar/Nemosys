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
#ifndef NEMOSYS_PATRAN_H_
#define NEMOSYS_PATRAN_H_

#include "nemosys_export.h"
#include "Mesh/meshBase.H"

namespace PATRAN {

class NEMOSYS_EXPORT patran
{
  public:
    patran(std::shared_ptr<meshBase> fullMesh,
           const std::string &inFnameVtk,
           const std::string &outFnameNeu,
           const std::map<int, int> &faceTypeMap,
           const std::map<int, int> &nodeTypeMap,
           const std::map<int, bool> &nodeStructuralMap,
           const std::map<int, bool> &nodeMeshMotionMap,
           const std::map<int, bool> &nodeThermalMap,
           const std::vector<int> &nppItr);

  public:
    void write25(std::ofstream &outputStream) const;
    void write26(std::ofstream &outputStream) const;
    void write1(std::ofstream &outputStream) const;
    void write2(std::ofstream &outputStream) const;
    void write6(std::ofstream &outputStream);
    void write8(std::ofstream &outputStream);
    void write99(std::ofstream &outputStream) const;
    bool comparePatch(int i, int j);

  private:
    std::shared_ptr<meshBase> fullMesh;
    std::string inFnameVtk;
    std::string outFnameNeu;
    std::shared_ptr<meshBase> volMeshBase;
    std::shared_ptr<meshBase> surfMeshBase;
    std::map<int, int> faceTypeMap;
    std::map<int, int> nodeTypeMap;
    std::map<int, bool> nodeStructuralMap;
    std::map<int, bool> nodeMeshMotionMap;
    std::map<int, bool> nodeThermalMap;
    std::vector<int> nppVec;
    std::map<int, std::string> face2nodes;
    std::map<nemId_t, std::vector<int>> boundaryNodeId2PatchNo;

};

} // namespace PATRAN

#endif  // NEMOSYS_PATRAN_H_
