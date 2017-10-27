#include <meshPhys.H>

// constructor for PointDataArray
// Moved from header to ensure access to flatten function
PointDataArray::PointDataArray( std::string _name, int _numComponent, int _numTuple, 
                    const std::vector<std::vector<double>>& _pntData):
                    name(_name), numComponent(_numComponent), numTuple(_numTuple)
                    { pntData = flatten (_pntData) ; };

// computes the gradient of point data at a cell using 
// derivatives of shape interpolation functions
std::vector<double> meshPhys::ComputeGradAtCell(int cell, int array)
{
  if (!pointData)
    pointData = dataSet->GetPointData();   
  if (pointData)
  {

    vtkIdList* point_ids = dataSet->GetCell(cell)->GetPointIds();
    int numPointsInCell = point_ids->GetNumberOfIds();
    int dim = pntData[array].getNumComponent();
    
    // populating array with point data of cell
    double values[dim*numPointsInCell];
    for (int i = 0; i < numPointsInCell; ++i)
    {
			int id = (int) point_ids->GetId(i);
      for (int j = 0; j < dim; ++j)
        values[i*dim +j] = pntData[array](id,j); 
    }
   
    double derivs[dim*3];
    double tmp[3];
    // getting gradient of field over cell (jacobian matrix for data)
    dataSet->GetCell(cell)->Derivatives(0,tmp,values,dim,derivs); 
    /* The Derivatives member for a cell computes the inverse of the jacobian
       transforming physical coordinates to parametric space, the derivatives of
       the shape functions in parametric space and the interpolated values of 
       the derivatives of data for the cell in parametric space. It then computes 
       the matrix product of the inverse jacobian with the interpolated derivatives 
       to transform them to physical coordinates */
    
    std::vector<double> gradient(derivs, derivs+dim*dim); 
    return gradient;
  }
}

// computes value of point data at a cell center using shape interpolation functions
std::vector<double> meshPhys::ComputeValAtCell(int cell, int array)
{
  if (!pointData)
    pointData = dataSet->GetPointData();
  if (pointData)
  {
    vtkIdList* point_ids = dataSet->GetCell(cell)->GetPointIds();
    int numPointsInCell = point_ids->GetNumberOfIds();
    int dim = pntData[array].getNumComponent();
    
    // evaluate center of cell in cartesian coords
    std::vector<double*> cell_coords = getCellCoords(cell,numPointsInCell);
    double x, y, z;
    x = y = z = 0;
    for (int i = 0; i < numPointsInCell; ++i)
    {
      x+= cell_coords[i][0];
      y+= cell_coords[i][1];
      z+= cell_coords[i][2];
      delete cell_coords[i];   
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
    dataSet->GetCell(cell)->EvaluatePosition(center, NULL, subId, pcoords, minDist2, weights);    
    /* The EvaluatePosition member for a cell takes the center in cartesian coordinates
       and evaluates its location in parametric space as well as the interpolation weights
       at that point. The NULL parameter prevents computation of closest point on cell
       and distance from center to cell (because center is in cell already) */

    // compute value of data at center of cell
    std::vector<double> values(dim); 
    for (int i = 0; i < dim; ++i) // loop over values per vertex
    {
      double val=0;
      for (int j = 0; j < numPointsInCell; ++j) // loop over vertices in cell
      {
        int id = (int) point_ids->GetId(j);
        val += pntData[array](id,i)*weights[j];
      }
      values[i] = val;
    }
    return values;
  }
}

// compute L2 norm of gradient of point data at each cell 
std::vector<double> meshPhys::ComputeL2GradAtAllCells(int array)
{
  std::vector<double> result;
  result.resize(numberOfCells);
  for (int i = 0; i < numberOfCells; ++i)
  {
    result[i] = L2_Norm(ComputeGradAtCell(i, array)); 
  }

  return result;
}

// compute value of point data at center of each cell
std::vector<double> meshPhys::ComputeValAtAllCells(int array)
{
  int dim = getDimArray(array);
  std::vector<double> result(numberOfCells*dim);
  for (int i = 0; i < numberOfCells; ++i)
  {
    std::vector<double> values = ComputeValAtCell(i, array);
    for (int j = 0; j < dim; ++j)
      result[i*dim + j] = values[j];
  }
  return result;
}

// compute L2 norm of value of point data at center of each cell
std::vector<double> meshPhys::ComputeL2ValAtAllCells(int array)
{
  std::vector<double> result(numberOfCells); 
  for (int i = 0; i < numberOfCells; ++i)
    result[i] = L2_Norm(ComputeValAtCell(i, array));
  return result;
}

// get diameter of circumsphere of each cell
std::vector<double> meshPhys::GetCellLengths()
{
  std::vector<double> result;
  result.resize(numberOfCells);
  for (int i = 0; i < numberOfCells; ++i)
  {
    result[i] = std::sqrt(dataSet->GetCell(i)->GetLength2());
  } 
  return result;
}

// writes a background mesh with sizes defined by 
// point data interpolated to cell center
void meshPhys::writeBackgroundMSH(string filename, std::vector<double> sizes)
{
  std::ofstream outputStream(filename.c_str());
  if(!outputStream.good()) {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
 
  if (!dataSet) {
    std::cout << "No data to write" << std::endl;
    exit(2);
  } 
  // ---------  writing gmsh header ----------- //
  outputStream << "$MeshFormat" << std::endl
               << "2.2 0 8" << std::endl
               << "$EndMeshFormat" << std::endl; 

  // ---- get number of points and number of elements ---- //
  getNumberOfCells();
  getNumberOfPoints();

  // -------- ensure all cell types are tri/tet or below -------------- //
  int num_bad = 0;
  for (int i = 0; i < numberOfCells; i++)
  {
    int type_id = dataSet->GetCellType(i);
    if (!(type_id == 3 || type_id == 5 || type_id == 10))
    {
      std::cout << "Error: Only tetrahedral and triangular" 
                << " meshes can be written to gmsh format" << std::endl;
      exit(3);
    }
    if (!(type_id == 10))
      num_bad+=1;
  }

  // ------------------------ write point coords -------------------------- //
  outputStream << "$Nodes" << std::endl << numberOfPoints << std::endl;
  for (int i = 0; i < numberOfPoints; ++i)
  {
    double* pntcrds = getPointCoords(i);
    outputStream << i + 1 << " "
                 << pntcrds[0] << " "
                 << pntcrds[1] << " "
                 << pntcrds[2] << " " << std::endl;
  }
  outputStream << "$EndNodes" << std::endl;

  // ------------- write element type and connectivity --------------------- //
  outputStream << "$Elements" << std::endl << numberOfCells-num_bad << std::endl;
  //int k = 0;
  for (int i = 0; i < numberOfCells; ++i)
  {

    vtkIdList* point_ids = dataSet->GetCell(i)->GetPointIds();
    int numComponent = point_ids->GetNumberOfIds();
    int type_id = dataSet->GetCellType(i);
    if (type_id == 10)
    {
      outputStream << i + 1 << " ";
      switch(numComponent)
      {
        case 2:
        {
          break;
        }
        case 3:
        {
          outputStream << 2 << " " << 2 << " " << 1 << " " << 1 << " ";
          break;
        }
        case 4:
        {
          outputStream << 4 << " " << 2 << " " << 1 << " " << 1 << " ";
          break;
        }
      
        default: 
        {  
          std::cerr << "Components in cell should be less than 4"<< std::endl;
          exit(1);
        }
      }
      for (int j = 0; j < numComponent; ++j)
         outputStream << point_ids->GetId(j) + 1 << " ";
      outputStream << std::endl;
      //k+=1;
    }
  }
  outputStream << "$EndElements" << std::endl;
  // -------------------------- write cell data ---------------------------- // 

  std::string tmpname = "BackgroundSF";
  outputStream << "$ElementData" << std::endl
               << 1 << std::endl // 1 string tag
               << "\"" << tmpname // name of view
               << "\"" << std::endl 
               << 0 << std::endl // 0 real tag
               << 3 << std::endl // 3 int tags (dt index, dim of field, number of fields)
               << 0 << std::endl // dt index
               << 1 << std::endl // dim of field
               << numberOfCells-num_bad << std::endl; // number of fields
  //int i = 0;
  for (int j = 0; j < numberOfCells; ++j)
  {
    int type_id = dataSet->GetCellType(j);
    if (type_id == 10) {
      outputStream << j+1 << " ";
      outputStream << sizes[j] << " ";
      outputStream << std::endl;
    //  i+=1;
    }
  }
  outputStream << "$EndElementData" << std::endl;    

    

}


// --------------------------- Auxilliary Functions -------------------------//

// flattens vector of vectors
template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& v)
{
	std::size_t size = 0;
	for (const auto& sub : v)
		size+= sub.size();
	std::vector<T> result;
	result.reserve(size);
	for (const auto& sub : v)
		result.insert(result.end(), sub.begin(), sub.end());
	return result;
}

// compute L2 norm of vec
double L2_Norm(const std::vector<double>& x)
{
  double result = 0.0;
  for (int i = 0; i < x.size(); ++i)
    result += x[i]*x[i];
  
  return std::sqrt(result);
}

// adds two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& x, 
                         const std::vector<T>& y)
{
  if (x.size() != y.size())
  {
    std::cout << "vectors must be same length for addition" << std::endl;
    exit(1);
  }
  
  std::vector<T> result;
  result.reserve(x.size());
  std::transform(x.begin(), x.end(), y.begin(), 
                 std::back_inserter(result), std::plus<T>());
  return result;
} 




/* computes gradient at point by averaging gradients at cells 
// surrounding point
std::vector<double> meshPhys::ComputeGradAtPoint(int pnt, int array)
{
  if (!pointData)
    pointData = dataSet->GetPointData();
  if (pointData)
  {
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  	vtkIdType pntId = pnt;
    dataSet->GetCellPoints(pntId, cellIds);
		std::cout << "number of id: " << cellIds->GetNumberOfIds() << std::endl;
    for (int i = 0; i < cellIds->GetNumberOfIds(); ++i)
			std::cout << cellIds->GetId(i) << " ";
		std::cout << std::endl;
	//  double* pntCoords = getPointCoords(pnt);
  //  std::cout << pntCoords[0] << " " << pntCoords[1] << " " << pntCoords[2] << std::endl; 
  //  std::cout << "cell number: " << dataSet->FindCell(pntCoords, NULL, NULL, 1e-10,
  
	}
 
  std::vector<double> tmp;
  return tmp; 

}*/
