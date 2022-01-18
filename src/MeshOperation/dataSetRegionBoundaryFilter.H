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
/*
 * Inspired by vtkDataSetRegionSurfaceFilter, but with an algorithm using STL
 * instead of vtkFastGeomQuadStruct.
 */

#ifndef NEMOSYS_DATASETREGIONBOUNDARYFILTER_H_
#define NEMOSYS_DATASETREGIONBOUNDARYFILTER_H_

#include "nemosys_export.h"
#include <vtkPolyDataAlgorithm.h>

#include <string>

namespace NEM {
namespace MSH {

/**
 * @brief Like vtkDataSetRegionSurfaceFilter, this filter extracts boundaries of
 * materials, including both interfaces and external boundaries. Unlike
 * vtkDataSetRegionSurfaceFilter, can extract boundary edges for 2d data sets.
 * @details Output edges/faces are represented once for interfaces, with the
 * orientation chosen by the input cell with lower material number. If input is
 * a vtkPointSet, output points are shallow copy; otherwise, input points are
 * deep copy - in either case, all points are copied, even those unused in
 * output. Output point data is shallow copy of input point data; input cell
 * data is not copied. Unlike vtkDataSetRegionSurfaceFilter, point ordering
 * should match the edge/face described in OrigCellArrayName and
 * CellSideArrayName arrays.
 */
class NEMOSYS_EXPORT dataSetRegionBoundaryFilter : public vtkPolyDataAlgorithm {
 public:
  static dataSetRegionBoundaryFilter *New();
  vtkTypeMacro(dataSetRegionBoundaryFilter, vtkPolyDataAlgorithm)

  //@{
  /**
   * @brief If 2-D, finds edges; if 3-D, finds faces. Default is 3.
   */
  vtkSetClampMacro(Dimension, int, 2, 3);
  vtkGetMacro(Dimension, int);
  //@}

  //@{
  /**
   * @brief Input material/region array name. Default is "Material".
   * @details Assumed to be vtkIntArray in input cell data. Interfaces between
   * materials are included in output. Should not use the value -1.
   */
  vtkSetMacro(MaterialArrayName, const std::string &);
  vtkGetMacro(MaterialArrayName, const std::string &);
  //@}

  //@{
  /**
   * @brief Output region array name. Default is "Region".
   * @details If not empty, array will be added as vtkIntArray in output cell
   * data. To recover relation to material, see RegionToMaterialArrayName
   */
  vtkSetMacro(BoundaryRegionArrayName, const std::string &);
  vtkGetMacro(BoundaryRegionArrayName, const std::string &);
  //@}

  //@{
  /**
   * @brief Output region to material array name. Default is "RegionToMaterial".
   * @details If not empty, array will be added as vtkIntArray with two
   * components in output field data. Indexed by region (value in
   * BoundaryRegionArrayName array); value corresponds to material(s), with -1
   * corresponding to no material (for regions bounding only one material).
   */
  vtkSetMacro(RegionToMaterialArrayName, const std::string &);
  vtkGetMacro(RegionToMaterialArrayName, const std::string &);
  //@}

  //@{
  /**
   * @brief Output original cell id array name. Default is "OrigCellIds".
   * @details If not empty, array will be added as vtkIdTypeArray with two
   * components in output cell data. For each output edge/face, contains the
   * index of the input cell. Second component is -1 if external edge/face.
   */
  vtkSetMacro(OrigCellArrayName, const std::string &);
  vtkGetMacro(OrigCellArrayName, const std::string &);
  //@}

  //@{
  /**
   * @brief Output original cell edge/face id array name. Default is
   * "CellFaceIds".
   * @details If not empty, array will be added as vtkIntArray with two
   * components in output cell data. For each output edge/face, contains the
   * index of the corresponding edge/face of the cell given by the
   * OrigCellArrayName array. Second component is -1 if external edge/face.
   */
  vtkSetMacro(CellSideArrayName, const std::string &);
  vtkGetMacro(CellSideArrayName, const std::string &);
  //@}

 protected:
  dataSetRegionBoundaryFilter() = default;

  int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  int FillInputPortInformation(int port, vtkInformation *info) override;

  int Dimension{3};
  std::string MaterialArrayName{"Material"};
  std::string RegionToMaterialArrayName{"RegionToMaterial"};
  std::string OrigCellArrayName{"OrigCellIds"};
  std::string CellSideArrayName{"CellFaceIds"};
  std::string BoundaryRegionArrayName{"Region"};
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_DATASETREGIONBOUNDARYFILTER_H_
