#include "SizeFieldGeneration/ValSizeField.H"

#include <vtkCellData.h>
#include <vtkCellIterator.h>
#include <vtkDoubleArray.h>
#include <vtkIdList.h>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace ADP {

// constructor
ValSizeField::ValSizeField(vtkDataSet *_ds, int arrayID, double _dev_mult,
                           bool _maxIsmin)
    : SizeFieldBase(_ds, arrayID, _dev_mult, _maxIsmin, "ValueSF") {
  std::cout << "ValSizeField constructed" << std::endl;
}

// computes value of point data at a cell center using average of data
// at points defining cell
// NOTE: averaging is equivalent to using interpolation weights for cell center
std::vector<double> ValSizeField::computeValAtCell(vtkIdList *cell_point_ids,
                                                   vtkDataArray *da) {
  int dim = da->GetNumberOfComponents();
  std::vector<double> values(dim, 0.0);
  int numPointsInCell = cell_point_ids->GetNumberOfIds();
  // compute value of data at center of cell
  for (int j = 0; j < numPointsInCell; ++j) {
    std::vector<double> comps(dim);
    da->GetTuple(cell_point_ids->GetId(j), comps.data());
    for (int i = 0; i < dim; ++i) {
      values[i] += comps[i] / numPointsInCell;
    }
  }

  return values;
}

// compute value of point data at center of each cell
std::vector<std::vector<double>> ValSizeField::computeValAtAllCells(
    vtkDataSet *ds, vtkDataArray *da) {
  std::vector<std::vector<double>> result;
  result.reserve(ds->GetNumberOfCells());

  vtkCellIterator *it = ds->NewCellIterator();
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell()) {
    result.emplace_back(computeValAtCell(it->GetPointIds(), da));
  }
  it->Delete();

  return result;
}

// compute 2 norm of value of point data at center of each cell
std::vector<double> ValSizeField::computeL2ValAtAllCells(vtkDataSet *ds,
                                                         vtkDataArray *da) {
  std::vector<double> result;
  result.reserve(ds->GetNumberOfCells());

  vtkCellIterator *it = ds->NewCellIterator();
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell()) {
    result.emplace_back(
        nemAux::l2_Norm(computeValAtCell(it->GetPointIds(), da)));
  }
  it->Delete();

  return result;
}

// compute size field and insert as cell data into mesh's dataSet
void ValSizeField::computeSizeField(vtkDataArray *da) {
  // populate vector with 2 norm of gradient/value of physical variable
  std::vector<double> values = computeL2ValAtAllCells(ds, da);

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
