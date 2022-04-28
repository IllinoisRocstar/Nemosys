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
#ifndef NEMOSYS_NUCMESHSRV_H_
#define NEMOSYS_NUCMESHSRV_H_

#include "nemosys_export.h"
#include "Services/srvBase.H"

#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <SMESH_Gen.hxx>

#include "NucMesh/ShapeBase.H"

namespace NEM {
namespace SRV {

struct NEMOSYS_EXPORT NucMeshConf {
  /**
   * @brief Geometry and mesh. Shapes in the list are successively overlaid on
   * previous shapes.
   * @details Using shared_ptr instead of unique_ptr because SWIG does not
   * deal well with move-only types
   */
  std::vector<std::shared_ptr<NEM::NUCMESH::ShapeBase>> geometryAndMesh;
  /**
   * @brief Extrude the final 2d mesh along the z-axis
   * @details Each entry in @c extrudeSteps denotes a layer with that height.
   * Empty array means that final mesh will be 2d (lying in xy plane)
   */
  std::vector<double> extrudeSteps;
  /**
   * @brief File name of output geometry file
   */
  std::string geoFileName;
  /**
   * @brief Maximum mesh edge length
   */
  double maxMeshSize = 1.e12;
  /**
   * @brief Minimum mesh edge length
   */
  double minMeshSize = 1.e12;
  /**
   * @brief Boolean for NETGEN engine
   */
  bool useNetgen = false;


  /**
   * Helper to insert shapes into @c geometryAndMesh
   * @tparam Shape type derived from NEM::NUCMESH::ShapeBase
   * @param shape object to copy or move
   * @return reference to inserted shape (note not the same object passed in)
   */
  template <typename Shape>
  typename std::add_lvalue_reference<Shape>::type insertShape(Shape &&shape) {
    geometryAndMesh.emplace_back(std::shared_ptr<NEM::NUCMESH::ShapeBase>{
        new typename std::decay<Shape>::type(std::forward<Shape>(shape))});
    return dynamic_cast<Shape &>(*geometryAndMesh.back());
  }

  /**
   * Helper to construct and insert shapes into @c geometryAndMesh
   * @tparam Shape type derived from NEM::NUCMESH::ShapeBase
   * @tparam Args types of a constructor of @p Shape
   * @param args arguments to pass to constructor of @p Shape
   * @return reference to inserted shape
   */
  template <typename Shape, typename... Args>
  Shape &makeShape(Args &&...args) {
    geometryAndMesh.emplace_back(std::shared_ptr<NEM::NUCMESH::ShapeBase>{
        new Shape(std::forward<Args>(args)...)});
    return dynamic_cast<Shape &>(*geometryAndMesh.back());
  }

  std::shared_ptr<SMESH_Gen> generator{new SMESH_Gen{}};
};

class NEMOSYS_EXPORT NucMeshSrv : public srvBase {
 public:
  static NucMeshSrv *New();
  void SetConfiguration(const NucMeshConf &configuration);
  const NucMeshConf &GetConfiguration() const;
  vtkTypeMacro(NucMeshSrv, srvBase)

 protected:
  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  NucMeshSrv();
  NucMeshConf conf_;
};

}  // namespace SRV
}  // namespace NEM

#endif  // NEMOSYS_NUCMESHSRV_H_
