#include "GradSizeField.H"

#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellIterator.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace ADP {

// constructor
GradSizeField::GradSizeField(vtkDataSet *_ds, int arrayID, double _dev_mult,
                             bool _maxIsmin)
    : SizeFieldBase(_ds, arrayID, _dev_mult, _maxIsmin, "GradientSF") {
  std::cout << "GradSizeField constructed" << std::endl;
}

// computes the gradient of point data at a cell using
// derivatives of shape interpolation functions
std::vector<double> GradSizeField::computeGradAtCell(vtkCell *cell,
                                                     vtkDataArray *da) {
  vtkIdList *point_ids = cell->GetPointIds();
  vtkIdType numPointsInCell = point_ids->GetNumberOfIds();
  int dim = da->GetNumberOfComponents();

  // populating array with point data of cell
  std::vector<double> values(dim * numPointsInCell);
  for (vtkIdType i = 0; i < numPointsInCell; ++i) {
    std::vector<double> comps(dim);
    da->GetTuple(point_ids->GetId(i), comps.data());
    for (int j = 0; j < dim; ++j) values[i * dim + j] = comps[j];
  }
  // # vals per vertex * # deriv directions (x,y,z)
  std::vector<double> gradient(dim * 3);
  std::vector<double> tmp(3, 0.0);
  // getting gradient of field over cell (Jacobian matrix for data)
  cell->Derivatives(0, tmp.data(), values.data(), dim, gradient.data());
  /* The Derivatives member for a cell computes the inverse of the Jacobian
   * transforming physical coordinates to parametric space, the derivatives of
   * the shape functions in parametric space and the interpolated values of
   * the derivatives of data for the cell in parametric space. It then
   * computes the matrix product of the inverse Jacobian with the interpolated
   * derivatives to transform them to physical coordinates. */
  return gradient;
}

// compute 2 norm of gradient of point data at each cell
std::vector<double> GradSizeField::computeL2GradAtAllCells(vtkDataSet *ds,
                                                           vtkDataArray *da) {
  std::vector<double> result;
  result.reserve(ds->GetNumberOfCells());

  vtkCellIterator *it = ds->NewCellIterator();
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell()) {
    it->GetCell(cell);
    result.emplace_back(nemAux::l2_Norm(computeGradAtCell(cell, da)));
  }
  it->Delete();

  return result;
}

// compute size field and insert as cell data into mesh's dataSet
void GradSizeField::computeSizeField(vtkDataArray *da) {
  // populate vector with 2 norm of gradient/value of physical variable
  std::vector<double> values = computeL2GradAtAllCells(ds, da);

  if (values.empty()) {
    std::cerr << "size array hasn't been populated!" << std::endl;
    exit(1);
  }

  mutateValues(values);

  vtkSmartPointer<vtkDoubleArray> da2 = vtkSmartPointer<vtkDoubleArray>::New();
  da2->SetName(sfname.c_str());
  da2->SetNumberOfComponents(1);
  for (vtkIdType i = 0; i < ds->GetNumberOfCells(); ++i)
    da2->InsertNextTypedTuple(&values[i]);
  ds->GetCellData()->AddArray(da2);
}

}  // namespace ADP
}  // namespace NEM
