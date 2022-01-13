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
#ifndef NEMOSYS_SMESHGEOMESH_H_
#define NEMOSYS_SMESHGEOMESH_H_

#include "nemosys_export.h"
#include "Mesh/geoMeshBase.H"

#include <memory>

class SMESH_Gen;
class SMESH_Mesh;

namespace NEM {
namespace MSH {

class NEMOSYS_EXPORT smeshGeoMesh : public geoMeshBase {
 public:
  smeshGeoMesh();
  ~smeshGeoMesh() override;
  static smeshGeoMesh *New();
  vtkTypeMacro(smeshGeoMesh, geoMeshBase)

  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

  void setSMeshMesh(std::unique_ptr<SMESH_Mesh> &&mesh,
                    std::shared_ptr<SMESH_Gen> gen);
  const SMESH_Mesh &getSMESHMesh() const;

 private:
  void resetNative() override;

  static GeoMesh SmeshToGM(SMESH_Mesh &mesh);
  static void GMToSMESH(const GeoMesh& geoMesh, SMESH_Mesh &outMesh);

  std::shared_ptr<SMESH_Gen> gen_;
  std::unique_ptr<SMESH_Mesh> mesh_;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_SMESHGEOMESH_H_
