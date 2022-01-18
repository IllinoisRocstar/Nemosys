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
#ifndef NEMOSYS_VTKGEOMESH_H_
#define NEMOSYS_VTKGEOMESH_H_

#include "Mesh/geoMeshBase.H"

namespace NEM {
namespace MSH {

/**
 * @class vtkGeoMesh
 * A concrete implementation of @c geoMeshBase representing a mesh in a @c
 * vtkUnstructuredGrid
 */
class NEMOSYS_EXPORT vtkGeoMesh : public geoMeshBase {
 public:
  static vtkGeoMesh *New();
  vtkTypeMacro(vtkGeoMesh, geoMeshBase)

 public:
  /**
   * Create a vtkGeoMesh from a vtkUnstructuredGrid file.
   * @param fileName name of file
   * @param phyGrpArrayName name of a vtkDataArray in the CellData defining the
   * physical group of each cell
   * @return new vtkGeoMesh from file
   */
  static vtkGeoMesh *Read(const std::string &fileName,
                          const std::string &phyGrpArrayName = std::string());

 public:
  /**
   * Construct a vtkGeoMesh with an empty mesh
   */
  vtkGeoMesh();

  /**
   * Construct a vtkGeoMesh from an existing vtkUnstructuredGrid
   * @param inUnstructuredGrid existing vtkUnstructuredGrid; it is used, not
   * copied
   * @param phyGrpArrayName name of a vtkDataArray in @p inUnstructuredGrid 's
   * CellData defining the physical group of each cell
   */
  explicit vtkGeoMesh(vtkUnstructuredGrid *inUnstructuredGrid,
                      const std::string &phyGrpArrayName = std::string());
  ~vtkGeoMesh() override;

 public:
  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

 public:
  /**
   * Copy the mesh
   * @param dest copy the vtkUnstructuredGrid into @p dest
   */
  void getVtkMesh(vtkUnstructuredGrid *dest);
  /**
   * Set the vtkGeoMesh's mesh to @p vtkMesh
   * @param vtkMesh existing vtkUnstructured grid; note the pointer is used,
   * the vtkUnstructuredGrid is not copied
   */
  void setVtkMesh(vtkUnstructuredGrid *vtkMesh);

 private:
  /**
   * Create a @c GeoMesh from a @c vtkUnstructuredGrid
   * @details Assumes gmsh has already been initialized.
   * @param vtkMesh mesh
   * @param phyGrpArrayName name of cell data array describing physical group of
   * each cell
   * @return @c GeoMesh
   */
  static GeoMesh vtk2GM(vtkUnstructuredGrid *vtkMesh,
                        const std::string &phyGrpArrayName = std::string());

  void resetNative() override;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_VTKGEOMESH_H_
