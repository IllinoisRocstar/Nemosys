#include <GradSizeField.H>

// constructor
GradSizeField::GradSizeField(meshBase* _mesh, int arrayID,double _dev_mult, bool _maxIsmin)
{
  // setting private vars
  mesh = _mesh;
  dev_mult = _dev_mult;
  maxIsmin = _maxIsmin;
  // checking for point data
  int numArr = mesh->getDataSet()->GetPointData()->GetNumberOfArrays();
  if (arrayID >= numArr)
  {
    std::cout << "ERROR: arrayID is out of bounds" << std::endl;
    std::cout << "There are " << numArr << " point data arrays" << std::endl;
    exit(1);
  }
  else if (numArr < 1)
  {
    std::cout << "no point data found" << std::endl;
    exit(1);
  }
  // setting data array member
  da = mesh->getDataSet()->GetPointData()->GetArray(arrayID); 
  // setting name of size field
  std::string array_name = mesh->getDataSet()->GetPointData()->GetArrayName(arrayID);
  int dim = da->GetNumberOfComponents(); 
  sfname = array_name.append("GradientSF");
  
  { // checking for name conflicts and removing SF with same name if it exists
    vtkCellData* cd = mesh->getDataSet()->GetCellData();
    if (cd->GetNumberOfArrays())
    { 
      for (int i = 0; i < cd->GetNumberOfArrays(); ++i)
      {
        std::string currname = cd->GetArrayName(i);
        if (!sfname.compare(currname))
        {
          std::cout << "Found size field identifier in cell data: " << currname << std::endl;
          std::cout << "Removing " << currname << " from dataSet" << std::endl;
          mesh->unsetCellDataArray(i);
          break;
        }
      }
    }
  }
  std::cout << "GradSizeField constructed" << std::endl; 
}

// computes the gradient of point data at a cell using 
// derivatives of shape interpolation functions
std::vector<double> GradSizeField::computeGradAtCell(int cell, int array)
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
    
    // populating array with point data of cell
    double values[dim*numPointsInCell];
    for (int i = 0; i < numPointsInCell; ++i)
    {
      int id = (int) point_ids->GetId(i);
      double comps[dim];
      da->GetTuple(id,comps);
      for (int j = 0; j < dim; ++j)
        values[i*dim +j] = comps[j]; 
    }
    double derivs[dim*3]; // # vals per vertex * # deriv directions (x,y,z)
    double tmp[3];
    // getting gradient of field over cell (jacobian matrix for data)
    mesh->getDataSet()->GetCell(cell)->Derivatives(0,tmp,values,dim,derivs); 
    /* The Derivatives member for a cell computes the inverse of the jacobian
       transforming physical coordinates to parametric space, the derivatives of
       the shape functions in parametric space and the interpolated values of 
       the derivatives of data for the cell in parametric space. It then computes 
       the matrix product of the inverse jacobian with the interpolated derivatives 
       to transform them to physical coordinates */
    
    std::vector<double> gradient(derivs, derivs+dim*3); 
    return gradient;
  
  }
  else
  {
    std::cout << "no point data found" << std::endl;
    exit(1);
  }
}

// compute 2 norm of gradient of point data at each cell
std::vector<double> GradSizeField::computeL2GradAtAllCells(int array)
{
  std::vector<double> result(mesh->getNumberOfCells()); 
  for (int i = 0; i < mesh->getNumberOfCells(); ++i)
  { 
    result[i] = l2_Norm(computeGradAtCell(i, array));
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
    std::cout << "size array hasn't been populated!" << std::endl;
    exit(1);
  }

  mutateValues(values);  
  mesh->setCellDataArray(&sfname[0u], values);
  mesh->setSFBool(1);
}
