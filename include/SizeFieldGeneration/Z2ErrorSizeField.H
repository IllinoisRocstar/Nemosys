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
#ifndef NEMOSYS_Z2ERRORSIZEFIELD_H_
#define NEMOSYS_Z2ERRORSIZEFIELD_H_

#include "nemosys_export.h"
#include "SizeFieldGeneration/SizeFieldGen.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT Z2ErrorSizeField : public SizeFieldBase {
  // constructors and destructors
 public:
  Z2ErrorSizeField(vtkDataSet *_ds, int arrayID, int _order);

  ~Z2ErrorSizeField() override {
    std::cout << "Z2ErrorSizeField destroyed" << std::endl;
  }

  // methods
 public:
  // perform patch recovery to get RMSE of recovered values of arrayID
  double computeNodalError(int arrayID) const;
  // compute size field and insert as cell data into mesh's dataSet
  void computeSizeField(vtkDataArray *da) override;

 public:
  int getOrder() const { return order; }
  void setOrder(int _order) { order = _order; }

 private:
  int order;
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_Z2ERRORSIZEFIELD_H_
