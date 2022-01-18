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
#ifndef NEMOSYS_VALSIZEFIELD_H_
#define NEMOSYS_VALSIZEFIELD_H_

#include "nemosys_export.h"
#include "SizeFieldGeneration/SizeFieldGen.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT ValSizeField : public SizeFieldBase {
  // constructors and destructors
 public:
  ValSizeField(vtkDataSet *_ds, int arrayID, double _dev_mult, bool _maxIsmin);

  ~ValSizeField() override {
    std::cout << "ValSizeField destroyed" << std::endl;
  }

  // methods
 public:
  // computes value of point data at a cell center using average of data
  // at points defining cell
  static std::vector<double> computeValAtCell(vtkIdList *cell_point_ids,
                                              vtkDataArray *da);
  // compute value of point data at center of each cell
  static std::vector<std::vector<double>> computeValAtAllCells(
      vtkDataSet *ds, vtkDataArray *da);
  // compute L2 norm of value of point data at center of each cell
  static std::vector<double> computeL2ValAtAllCells(vtkDataSet *ds,
                                                    vtkDataArray *da);
  // compute size field and insert as cell data into mesh's dataSet
  void computeSizeField(vtkDataArray *da) override;
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_VALSIZEFIELD_H_
