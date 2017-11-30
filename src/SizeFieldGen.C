#include <SizeFieldGen.H>

SizeFieldGen* SizeFieldGen::Create(meshBase* _mesh, std::string method, int arrayID,
                                   double _dev_mult, bool _maxIsmin)
{
  if (!method.compare("value"))
  {
    ValSF* valsf = new ValSF(_mesh, arrayID, _dev_mult, _maxIsmin);
    valsf->computeSizeField(arrayID);
    return valsf;
  }
  else if (!method.compare("gradient"))
  {
    GradSF* gradsf = new GradSF(_mesh, arrayID, _dev_mult, _maxIsmin);
    gradsf->computeSizeField(arrayID);
    return gradsf;
  }
  else if (!method.compare("error")){}
  else
  {
    std::cout << "Specified method " << method << " is not supported" << std::endl;
    std::cout << "Available methods are gradient, value and error" << std::endl;
    exit(1);
  }
  
}

GradSF::GradSF(meshBase* _mesh, int arrayID,double _dev_mult, bool _maxIsmin)
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
          mesh->unsetCellDataArray(i);
          break;
        }
      }
    }
  }
  std::cout << "GradSF constructed" << std::endl; 
}


ValSF::ValSF(meshBase* _mesh, int arrayID, double _dev_mult, bool _maxIsmin)
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
  sfname = array_name.append("ValueSF");
  
  { // checking for name conflicts and removing SF with same name if it exists
    vtkCellData* cd = mesh->getDataSet()->GetCellData();
    if (cd->GetNumberOfArrays())
    { 
      for (int i = 0; i < cd->GetNumberOfArrays(); ++i)
      {
        std::string currname = cd->GetArrayName(i);
        if (!sfname.compare(currname))
        {
          mesh->unsetCellDataArray(i);
          break;
        }
      }
    }
  }
  std::cout << "ValSF constructed" << std::endl; 
}

// computes the gradient of point data at a cell using 
// derivatives of shape interpolation functions
std::vector<double> GradSF::computeGradAtCell(int cell, int array)
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

 
// computes value of point data at a cell center using shape interpolation functions
std::vector<double> ValSF::computeValAtCell(int cell, int array)
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
    // evaluate center of cell in cartesian coords
    std::map<int, std::vector<double>> cell_coords = mesh->getCell(cell);
    std::map<int, std::vector<double>>::iterator it = cell_coords.begin();

    double x, y, z;
    x = y = z = 0;
    while (it != cell_coords.end())
    {
      x += it->second[0];
      y += it->second[1];
      z += it->second[2];
      ++it;
    }
    
    double center[3];
    center[0] = x/numPointsInCell;
    center[1] = y/numPointsInCell;
    center[2] = z/numPointsInCell;
    
    // evaluate parametric position of center and compute interpolation weights
    int subId;
    double minDist2; // not used
    double pcoords[3];
    double weights[numPointsInCell];
    mesh->getDataSet()->GetCell(cell)->
      EvaluatePosition(center, NULL, subId, pcoords, minDist2, weights);    
    /* The EvaluatePosition member for a cell takes the center in cartesian coordinates
       and evaluates its location in parametric space as well as the interpolation weights
       at that point. The NULL parameter prevents computation of closest point on cell
       and distance from center to cell (because center is in cell already) */

    // compute value of data at center of cell
    std::vector<double> values(dim,0); 
/*    for (int i = 0; i < dim; ++i) // loop over values per vertex
    {
      double val=0;
      for (int j = 0; j < numPointsInCell; ++j) // loop over vertices in cell
      {
        int id = (int) point_ids->GetId(j);
        double comps[dim];
        da->GetTuple(id, comps);
        val += comps[i]*weights[j];
      }
      values[i] = val;
    }
*/
    for (int j = 0; j < numPointsInCell; ++j)
    {
      int id = point_ids->GetId(j);
      double comps[dim];
      da->GetTuple(id,comps);
      for (int i = 0; i < dim; ++i)
      {
        values[i] += comps[i]*weights[j];
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
std::vector<std::vector<double>> ValSF::computeValAtAllCells(int arrayID)
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

// compute L2 norm of value of point data at center of each cell
std::vector<double> ValSF::computeL2ValAtAllCells(int array)
{
  std::vector<double> result(mesh->getNumberOfCells()); 
  for (int i = 0; i < mesh->getNumberOfCells(); ++i)
    result[i] = L2_Norm(computeValAtCell(i, array));
  return result;
}

// compute size field and insert as cell data into mesh's dataSet
void ValSF::computeSizeField(int arrayID)
{
  // populate vector with L2 norm of gradient/value of physical variable
  std::vector<double> values = computeL2ValAtAllCells(arrayID); 
  
  if (values.empty())
  {
    std::cout << "size array hasn't been populated!" << std::endl;
    exit(1);
  }

  mutateValues(values);  
  mesh->setCellDataArray(&sfname[0u], values);
}

// identifies cells to refine and mutates current size values
// into a compatible size field for the mesh
void SizeFieldGen::mutateValues(std::vector<double>& values)
{
  // get circumsphere diameter of all cells 
  std::vector<double> lengths = mesh->getCellLengths();
  // find minmax of diameters
  std::vector<double> lengthminmax = getMinMax(lengths);
  // redefine minmax values for appropriate size definition reference
  if (maxIsmin)
    lengthminmax[1] = lengthminmax[0];
  else
    lengthminmax[1] *= 0.65;
  lengthminmax[0] -= lengthminmax[0]/2.; 



  // get mean and stdev of values 
  std::vector<double> meanStdev = getMeanStdev(values);
  // get bool array of which cells to refine based on multiplier of stdev
  std::vector<int> cells2Refine = cellsToRefine(values, meanStdev[0]+meanStdev[1]*dev_mult);
  // normalize values by mean
  std::vector<double> values_norm = (1./meanStdev[0])*values;  
  // take the reciprocal of values for size definition (high value -> smaller size)
  reciprocal_vec(values);
  // scale values to min max circumsphere diam of cells 
  // now, values represents a size field
  std::vector<double> valuesMinMax = getMinMax(values);
  scale_vec_to_range(values, valuesMinMax, lengthminmax);
  // setting sizes (values) to f*max element diam based on return of cellsToRefine function
  for (int i = 0; i < values.size(); ++i)
  {
    if (!cells2Refine[i])
      values[i] = lengthminmax[1]; // if cell shouldn't be refined, size set to min of diams
  }
}
