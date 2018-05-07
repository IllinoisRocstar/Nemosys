#include <ValSizeField.H>

// constructor
ValSizeField::ValSizeField(meshBase* _mesh, int arrayID, double _dev_mult, bool _maxIsmin)
{
  initialize(_mesh, arrayID, _dev_mult, _maxIsmin, "ValueSF");
  std::cout << "ValSizeField constructed" << std::endl; 
}

// computes value of point data at a cell center using average of data
// at points defining cell
// NOTE: averaging is equivalent to using interpolation weights for cell center
std::vector<double> ValSizeField::computeValAtCell(int cell, int array)
{
  if (!mesh)
  {
    std::cout << "no mesh object exists!" << std::endl;
    exit(1);
  }
  
  if (da)
  {
    vtkIdList* point_ids = mesh->getDataSet()->GetCell(cell)->GetPointIds();
    int numPointsInCell = point_ids->GetNumberOfIds();
    int dim = da->GetNumberOfComponents();
    // compute value of data at center of cell
    std::vector<double> values(dim,0); 
    for (int j = 0; j < numPointsInCell; ++j)
    {
      int id = point_ids->GetId(j);
      double comps[dim];
      da->GetTuple(id,comps);
      for (int i = 0; i < dim; ++i)
      {
        values[i] += comps[i]/numPointsInCell;
      }
    }
    return values;
  }
  else
  {
    std::cout << "no point data found" << std::endl;
    exit(1);
  }
}

// compute value of point data at center of each cell
std::vector<std::vector<double>> ValSizeField::computeValAtAllCells(int arrayID)
{
  int dim = mesh->getDataSet()->GetPointData()->GetArray(arrayID)->GetNumberOfComponents();
  std::vector<std::vector<double>> result(mesh->getNumberOfCells());
  for (int i = 0; i < mesh->getNumberOfCells(); ++i)
  {
    result[i].resize(dim);
    result[i] = computeValAtCell(i, arrayID);
  }
  return result;
}

// compute 2 norm of value of point data at center of each cell
std::vector<double> ValSizeField::computeL2ValAtAllCells(int array)
{
  std::vector<double> result(mesh->getNumberOfCells()); 
  for (int i = 0; i < mesh->getNumberOfCells(); ++i)
    result[i] = l2_Norm(computeValAtCell(i, array));
  return result;
}

// compute size field and insert as cell data into mesh's dataSet
void ValSizeField::computeSizeField(int arrayID)
{
  // populate vector with 2 norm of gradient/value of physical variable
  std::vector<double> values = computeL2ValAtAllCells(arrayID); 
 
  if (values.empty())
  {
    std::cout << "size array hasn't been populated!" << std::endl;
    exit(1);
  }

  mutateValues(values);  
  mesh->setCellDataArray(&sfname[0u], values);
  mesh->setSFBool(1);
}
