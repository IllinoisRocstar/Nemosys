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
#ifndef NEMOSYS_REFINE_H_
#define NEMOSYS_REFINE_H_

#include <string>

#include <MAdLib.h>

#include "nemosys_export.h"
#include "Mesh/meshBase.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT Refine {
  // constructors and destructors
 public:
  Refine(meshBase *_mesh, const std::string &method, int arrayID,
         double dev_mult, bool maxIsmin, double edge_scale,
         const std::string &_ofname, double sizeFactor = 1.0, int order = 1);
  ~Refine();

  // methods
 public:
  void run(bool transferData, bool bndryConstraint = false);

 private:
  meshBase *mesh;      // mesh to be refined
  MAd::pMesh MadMesh;  // MAdLib mesh object generated from converted meshBase
  MAd::BackgroundSF *bSF;     // background sizeField
  MAd::PWLSField *pwlSF;      // piecewise linear size field
  MAd::MeshAdapter *adapter;  // adapter
  std::string ofname;
  MAd::pGEntity bnd;
  bool bndrConst;

  // helpers
 private:
  void initUniform(double edge_scale);
  void initAdaptive(int arrayID, const std::string &method);

  void classifyBoundaries();
  void unClassifyBoundaries();
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_REFINE_H_
