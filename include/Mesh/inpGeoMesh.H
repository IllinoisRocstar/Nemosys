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
#ifndef NEMOSYS_INPGEOMESH_H_
#define NEMOSYS_INPGEOMESH_H_

#include "Mesh/geoMeshBase.H"

#include <map>
#include <set>
#include <utility>

namespace NEM {
namespace MSH {

/**
 * @brief Class representing meshes in CalculiX input deck (similar to ABAQUS)
 * @details Currently supports the keywords ELEMENT ("C3D4", "C3D8", and "C3D6"
 * types, ELSET parameter), ELSET, NSET, and SURFACE (TYPE=ELEMENT only). By
 * default, physical groups with dimension equal to mesh dimension are turned
 * into ELSETs, and NSETs and physical groups with dimension one less than mesh
 * dimension are turned into SURFACEs and NSETs.
 */
class NEMOSYS_EXPORT inpGeoMesh : public geoMeshBase {
 private:
  struct NEMOSYS_NO_EXPORT InpSets {
    /**
     * @brief NSET keyword; indexing matches GeoMesh::mesh
     */
    std::map<std::string, std::set<vtkIdType>> nodeSets;
    /**
     * @brief ELSET keyword; indexing matches GeoMesh::mesh
     */
    std::map<std::string, std::set<vtkIdType>> elSets;
    /**
     * @brief SURFACE keyword; indexing matches GeoMesh::sideSet
     */
    std::map<std::string, std::set<vtkIdType>> surfaces;
  };

 public:
  static inpGeoMesh *New();
  static inpGeoMesh *Read(const std::string &fileName);
  vtkTypeMacro(inpGeoMesh, geoMeshBase)

  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

  decltype(InpSets::nodeSets) &getNSets() { return inpSets_.nodeSets; }
  decltype(InpSets::elSets) &getElSets() { return inpSets_.elSets; }
  decltype(InpSets::surfaces) &getSurfaces() { return inpSets_.surfaces; }

 protected:
  inpGeoMesh();
  explicit inpGeoMesh(const std::string &fileName);

  void reconstructGeo() override;

 private:
  /**
   * @brief Ctor meant for ease of delegating to geoMeshBase ctor
   * @param mesh from parsing a .inp mesh file
   */
  explicit inpGeoMesh(std::pair<GeoMesh, InpSets> mesh);
  void resetNative() override;
  static std::pair<GeoMesh, InpSets> inp2GM(const std::string &fileName);

  /**
   * @brief Holds data specific to inp format
   */
  InpSets inpSets_;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_INPGEOMESH_H_
