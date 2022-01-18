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
#include "SizeFieldGeneration/Z2ErrorSizeField.H"

#include <vtkCellData.h>
#include <vtkCellIterator.h>

#include "PatchRecovery/patchRecovery.H"

namespace NEM {
namespace ADP {

Z2ErrorSizeField::Z2ErrorSizeField(vtkDataSet *_ds, int arrayID, int _order)
    : SizeFieldBase(_ds, arrayID, 0.0, false, "Z2ErrorSF"), order(_order) {
  std::cout << "Z2ErrorSizeField constructed" << std::endl;
}

double Z2ErrorSizeField::computeNodalError(int arrayID) const {
  std::vector<int> arrayIDs(1);
  arrayIDs[0] = arrayID;
  std::unique_ptr<PatchRecovery> recoverObj =
      std::unique_ptr<PatchRecovery>(new PatchRecovery(ds, order, arrayIDs));
  return recoverObj->computeNodalError()[0][0];
}

void Z2ErrorSizeField::computeSizeField(vtkDataArray *da) {
  int arrayID;
  ds->GetPointData()->GetArray(da->GetName(), arrayID);
  double aveError = computeNodalError(arrayID) / ds->GetNumberOfCells();

  std::string errorName = da->GetName();
  errorName += "ErrorIntegral";

  vtkSmartPointer<vtkDataArray> errorIntegrals =
      ds->GetCellData()->GetArray(errorName.c_str());
  vtkSmartPointer<vtkDataArray> nodeSizes =
      ds->GetPointData()->GetArray("nodeSizes");

  std::vector<double> sizeField;
  sizeField.reserve(ds->GetNumberOfCells());

  vtkCellIterator *it = ds->NewCellIterator();
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell()) {
    double interpSize = 0.0;
    it->GetCell(cell);

    double center[3];
    auto *weights = new double[cell->GetNumberOfPoints()];
    int subId;    // dummy
    double x[3];  // dummy
    cell->GetParametricCenter(center);
    cell->EvaluateLocation(subId, center, x, weights);

    for (int j = 0; j < cell->GetNumberOfPoints(); ++j) {
      int pntID = cell->GetPointId(j);
      double size[1];
      nodeSizes->GetTuple(pntID, size);
      interpSize += size[0] * weights[j];
    }
    delete[] weights;

    double error[1];
    errorIntegrals->GetTuple(it->GetCellId(), error);
    sizeField.emplace_back(
        interpSize / (std::pow(error[0] / aveError, 1.0 / order)) * sizeFactor);
  }
  it->Delete();

  vtkSmartPointer<vtkDoubleArray> da2 = vtkSmartPointer<vtkDoubleArray>::New();
  da2->SetName(sfname.c_str());
  da2->SetNumberOfComponents(1);
  for (vtkIdType i = 0; i < ds->GetNumberOfCells(); ++i)
    da2->InsertNextTypedTuple(&sizeField[i]);
  ds->GetCellData()->AddArray(da2);
}

}  // namespace ADP
}  // namespace NEM
