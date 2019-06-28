#include "GradSizeField.H"
#include "AuxiliaryFunctions.H"

#include <vtkCell.h>

// constructor
GradSizeField::GradSizeField(meshBase *_mesh, int arrayID, double _dev_mult,
                             bool _maxIsmin)
{
  initialize(_mesh, arrayID, _dev_mult, _maxIsmin, "GradientSF");
  std::cout << "GradSizeField constructed" << std::endl;
}

// computes the gradient of point data at a cell using 
// derivatives of shape interpolation functions
std::vector<double> GradSizeField::computeGradAtCell(int cell, int array) const
{
  if (!mesh)
  {
    std::cerr << "no mesh object exists!" << std::endl;
    exit(1);
  }

  if (da)
  {
    vtkIdList *point_ids = mesh->getDataSet()->GetCell(cell)->GetPointIds();
    int numPointsInCell = point_ids->GetNumberOfIds();
    int dim = da->GetNumberOfComponents();

    // populating array with point data of cell
    auto *values = new double[dim * numPointsInCell];
    for (int i = 0; i < numPointsInCell; ++i)
    {
      int id = (int) point_ids->GetId(i);
      auto *comps = new double[dim];
      da->GetTuple(id, comps);
      for (int j = 0; j < dim; ++j)
        values[i * dim + j] = comps[j];
      delete[] comps;
    }
    // # vals per vertex * # deriv directions (x,y,z)
    auto derivs = new double[dim * 3];
    double tmp[3];
    // getting gradient of field over cell (Jacobian matrix for data)
    mesh->getDataSet()->GetCell(cell)->Derivatives(0, tmp, values, dim, derivs);
    /* The Derivatives member for a cell computes the inverse of the Jacobian
     * transforming physical coordinates to parametric space, the derivatives of
     * the shape functions in parametric space and the interpolated values of
     * the derivatives of data for the cell in parametric space. It then
     * computes the matrix product of the inverse Jacobian with the interpolated
     * derivatives to transform them to physical coordinates. */

    std::vector<double> gradient(derivs, derivs + dim * 3);
    delete derivs;
    delete[] values;
    return gradient;
  }
  else
  {
    std::cerr << "no point data found" << std::endl;
    exit(1);
  }
}

// compute 2 norm of gradient of point data at each cell
std::vector<double> GradSizeField::computeL2GradAtAllCells(int array) const
{
  std::vector<double> result(mesh->getNumberOfCells());
  for (int i = 0; i < mesh->getNumberOfCells(); ++i)
  {
    result[i] = nemAux::l2_Norm(computeGradAtCell(i, array));
  }
  return result;
}

// compute size field and insert as cell data into mesh's dataSet
void GradSizeField::computeSizeField(int arrayID)
{
  // populate vector with 2 norm of gradient/value of physical variable
  std::vector<double> values = computeL2GradAtAllCells(arrayID);

  if (values.empty())
  {
    std::cerr << "size array hasn't been populated!" << std::endl;
    exit(1);
  }

  mutateValues(values);
  mesh->setCellDataArray(sfname.c_str(), values);
  mesh->setSFBool(true);
}
