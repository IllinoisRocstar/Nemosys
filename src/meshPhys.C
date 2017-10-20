#include <meshPhys.H>

// compute gradient at center of cell
// for generality, we need all information of data array
vector<double> meshPhys::ComputeGradAtCell(int cell, int array)
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
   
    double derivs[dim*dim];
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

// computes gradient at point by averaging gradients at cells 
// surrounding point
std::vector<double> meshPhys::ComputeGradAtPoint(int pnt, int array)
{
  if (!pointData)
    pointData = dataSet->GetPointData();
  if (pointData)
  {
    vtkIdList* cellids = vtkIdList::New();
    dataSet->GetPointCells(pnt, cellids);
    std::cout << "number of id: " << cellids->GetId(0) << std::endl;
    /*double* pntCoords = getPointCoords(pnt);
    std::cout << pntCoords[0] << " " << pntCoords[1] << " " << pntCoords[2] << std::endl; 
    std::cout << "cell number: " << dataSet->FindCell(pntCoords, NULL, NULL, 1e-10,*/
  
  }
 
  std::vector<double> tmp;
  return tmp; 

}


// writes a background mesh with sizes as point data
void meshPhys::writeBackgroundMSH(string filename, std::vector<double> sizes)
{
    
  std::ofstream outputStream(filename.c_str());
  if(!outputStream.good()) 
  {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
  
  if (!dataSet) 
  {
    std::cout << "No data to write" << std::endl;
    exit(2);
  }
  
  // ---- get number of points and number of elements ---- //
  getNumberOfCells();
  getNumberOfPoints();
  
  if (sizes.size() != numberOfCells)
  {
    std::cout << "number of sizes is unequal to number of cells!" << std::endl; 
    exit(3);
  }
  

  // ---------  writing gmsh header ----------- //
  outputStream << "$MeshFormat" << std::endl
               << "2.2 0 8" << std::endl
               << "$EndMeshFormat" << std::endl; 

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
    if (!(type_id == 5 || type_id == 10))
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
  int k = 0;
  for (int i = 0; i < numberOfCells; ++i)
  {

    vtkIdList* point_ids = dataSet->GetCell(i)->GetPointIds();
    int numComponent = point_ids->GetNumberOfIds();
    int type_id = dataSet->GetCellType(i);
    if (type_id==5 || type_id == 10)
    {
      outputStream << k + 1 << " ";
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
      k+=1;
    }
  }
  outputStream << "$EndElements" << std::endl;
  
  // -------------------------- write cell data ---------------------------- // 

  std::string tmpname = "BackgroundSF";
  outputStream << "$ElementData" << std::endl
               << 1 << std::endl // 1 string tag
               << "\"" << (tmpname) // name of view
               << "\"" << std::endl 
               << 0 << std::endl // 0 real tag
               << 3 << std::endl // 3 int tags (dt index, dim of field, number of fields)
               << 0 << std::endl // dt index
               << 1 << std::endl // dim of field
               << numberOfCells-num_bad << std::endl; // number of fields
  int i = 0;
  for (int j = 0; j < numberOfCells; ++j)
  {
    int type_id = dataSet->GetCellType(j);
    if (type_id == 5 || type_id == 10) {
      outputStream << i+1 << " ";
      outputStream << sizes[j] << " ";
    outputStream << std::endl;
      i+=1;
    }
  }
  outputStream << "$EndElementData" << std::endl;    
}



 
        /*std::map<std::string, std::vector<std::vector<double>>>::iterator it = pntData.begin();
        while( it != pntData.end())
        {
          std::cout << it->first << ": " << (it->second)[1000][2] << std::endl;
          it++;
        }
    */
