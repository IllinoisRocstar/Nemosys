#include <meshPhys.H>

  /*  Algorithm: ComputeGradAtCenter:
      BEGIN:
        for each CELL in CELLS:
          JacboianInverse(double** inverse, double derivs[12]);

          |Now, we have the inverse jacobian (derivatives of shape function in iso-param space)
          |We use the interpolation derivatives (in iso-param space) to get the value of 
          |the functions derivative at cell center (in iso-param space)
          |Then, we use the inverse jacobian to transform back to the original coordinates

          |derivs[12] -> derivs[4][3] (deriv of shape function for each point)
          |Consider point data U=(Ux,Uy,Uz)
          
          dUx/dx = dUx/dy = dUx/dz = 0;
          for each POINT in CELL:
            dUx/dx += derivs[POINT][0]*ux[POINT]; 
            dUx/dy += derivs[POINT][1]*ux[POINT];
            dUx/dy += derivs[POINT][2]*ux[POINT]; 

          [dU/dx dU/dy dU/dz]_orig = inverse*[dU/dx dU/dy dU/dz];              
          return [dU/dx dU/dy dU/dz]_orig; 
      END
  */

// compute gradient at center of cell
// for generality, we need all information of data array
std::vector<std::vector<double>> meshPhys::ComputeGradAtCenter(int cell, int array)
{
  if (!pointData)
    pointData = dataSet->GetPointData();   
  if (pointData)
  {

    vtkIdList* point_ids = dataSet->GetCell(cellId)->GetPointIds();
    int numPointsInCell = point_ids->GetNumberOfIds();
    
    // computing inverse jacobian for element
    double** invJacobian;
    double  derivs[numPointsInCell*3];
    dataSet->GetCell(i)->JacobianInverse(invJacobian, derivs);
    
    double df = 0; 
    for (int i = 0; i < numPointsInCell; ++i)
    {
      vtkIdType id = point_ids->GetId(i);
      df += derivs 
      
    }
  
  } 

  }


}

// returns isoparametric coordinate (for generality)
// for 1st order tets it's always {0 0 0 1 0 0 0 1 0 0 0 1}
vector<double*> meshPhys::getParametricCoords()
{
  vector<double*> All_parametric_coords(numberOfCells);
  for (int i = 0; i < numberOfCells; ++i)
    All_parametric_coords[i] = dataSet->GetCell(i)->GetParametricCoords();
  return All_parametric_coords;
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
