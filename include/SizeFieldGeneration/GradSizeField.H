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
#ifndef NEMOSYS_GRADSIZEFIELD_H_
#define NEMOSYS_GRADSIZEFIELD_H_

#include "nemosys_export.h"
#include "SizeFieldGeneration/SizeFieldGen.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT GradSizeField : public SizeFieldBase {
  // constructors and destructors
 public:
  GradSizeField(vtkDataSet *_ds, int arrayID, double _dev_mult, bool _maxIsmin);

  ~GradSizeField() override {
    std::cout << "GradSizeField destroyed" << std::endl;
  }

  // methods
 public:
  // computes the gradient of point data at a cell using
  // derivatives of shape interpolation functions
  static std::vector<double> computeGradAtCell(vtkCell *cell, vtkDataArray *da);
  // compute L2 norm of gradient of point data at each cell
  static std::vector<double> computeL2GradAtAllCells(vtkDataSet *ds,
                                                     vtkDataArray *da);
  // compute size field and insert as cell data into mesh's dataSet
  void computeSizeField(vtkDataArray *da) override;
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_GRADSIZEFIELD_H_
