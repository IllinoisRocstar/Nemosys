#include <meshPhys.H>
// constructor for PointDataArray
// Moved from header to ensure access to flatten function
PointDataArray::PointDataArray( std::string _name, int _numComponent, int _numTuple, 
                    const std::vector<std::vector<double>>& _pntData):
                    name(_name), numComponent(_numComponent), numTuple(_numTuple)
                    { pntData = flatten (_pntData) ; }

std::vector<std::vector<double>> PointDataArray::getFoldData() 
{ 
  return fold(pntData, numComponent); 
}
  
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
   
    double derivs[dim*3]; // # vals per vertex * # deriv directions (x,y,z)
    double tmp[3];
    // getting gradient of field over cell (jacobian matrix for data)
    dataSet->GetCell(cell)->Derivatives(0,tmp,values,dim,derivs); 
    /* The Derivatives member for a cell computes the inverse of the jacobian
       transforming physical coordinates to parametric space, the derivatives of
       the shape functions in parametric space and the interpolated values of 
       the derivatives of data for the cell in parametric space. It then computes 
       the matrix product of the inverse jacobian with the interpolated derivatives 
       to transform them to physical coordinates */
    
    std::vector<double> gradient(derivs, derivs+dim*3); 
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

// generate background size field based on values or gradient
// and write it to msh file
void meshPhys::createSizeField(int array_id, std::string method, 
                               double dev_mult, bool maxIsmin)
{
  // get name of array for proper data naming in output
  std::string array_name =  getPointData(array_id).getName(); 
  // get dim of data (number of values per vertex)
  int dim = getDimArray(array_id);
  // get circumsphere diameter of all cells 
  std::vector<double> lengths = GetCellLengths();
  // find minmax of diameters
  std::vector<double> lengthminmax = getMinMax(lengths);
  // redefine minmax values for appropriate size definition reference
  if (maxIsmin)
    lengthminmax[1] = lengthminmax[0];
  else
    lengthminmax[1] *= 0.65;
  lengthminmax[0] -= lengthminmax[0]/2.; 

  // populate vector with L2 norm of gradient/value of physical variable
  std::vector<double> values;
 
  if (method.compare("grad") == 0)
    values = ComputeL2GradAtAllCells(array_id);
  else if (method.compare("val") == 0)
    values = ComputeL2ValAtAllCells(array_id); 
  else
  {
    std::cout << "Error creating size field" << std::endl
              << "Method must be \"val\" or \"grad\" " << std::endl;
    exit(1);
  }
  
  if (values.empty())
  {
    std::cout << "size array hasn't been populated!" << std::endl;
    exit(1);
  }

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

  // write cells2Refine to file
  std::ofstream outputstream("cells2Refine.txt");
  if (!outputstream.good())
  { 
    std::cout << "error opening file cells2Refine.txt" << std::endl;
  }
  for (int i = 0; i < cells2Refine.size(); ++i)
  {
    if (!cells2Refine[i])
    {
      vtkIdList* point_ids = dataSet->GetCell(i)->GetPointIds();
      outputstream << i << " ";
      for (int j = 0; j < point_ids->GetNumberOfIds(); ++j)
      {
        outputstream << point_ids->GetId(j) << " "; 
      }
      outputstream << std::endl;
    }
  }  
  // writing background mesh and taking care of non tet/tri entities
  writeBackgroundMSH("backgroundSF.msh", values);
}

// write cells to refine and connectivities to file
void meshPhys::writeCellsToRefine(int array_id, std::string method, double dev_mult)
{
  // get name of array for proper data naming in output
  std::string array_name =  getPointData(array_id).getName(); 

  // populate vector with L2 norm of gradient/value of physical variable
  std::vector<double> values;
 
  if (method.compare("grad") == 0)
    values = ComputeL2GradAtAllCells(array_id);
  else if (method.compare("val") == 0)
    values = ComputeL2ValAtAllCells(array_id); 
  else
  {
    std::cout << "Method must be \"val\" or \"grad\" " << std::endl;
    exit(1);
  }
  
  if (values.empty())
  {
    std::cout << "size array hasn't been populated!" << std::endl;
    exit(1);
  }

  // get mean and stdev of values 
  std::vector<double> meanStdev = getMeanStdev(values);
  // get bool array of which cells to refine based on multiplier of stdev
  std::vector<int> cells2Refine = cellsToRefine(values, meanStdev[0]+meanStdev[1]*dev_mult);

  // write cells2Refine to file
  std::ofstream outputstream("cells2Refine.txt");
  if (!outputstream.good())
  { 
    std::cout << "error opening file cells2Refine.txt" << std::endl;
  }
  for (int i = 0; i < cells2Refine.size(); ++i)
  {
    if (!cells2Refine[i])
    {
      vtkIdList* point_ids = dataSet->GetCell(i)->GetPointIds();
      outputstream << i << " ";
      for (int j = 0; j < point_ids->GetNumberOfIds(); ++j)
      {
        outputstream << point_ids->GetId(j) << " "; 
      }
      outputstream << std::endl;
    }
  }  
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
      
/* This function essentially does the following:

1) Registers the relevant data with the nodal data manager
2) Runs the adapter
3) Unclassifies the boundary elements for proper output
4) Writes the refined mesh to a .msh file without data
5) Gets the data after refinement
6) Converts the refined mesh to a vtk file without data
7) Writes the post-refinement data to a refined vtk and gmsh mesh 
8) Misc: Memory management, string trimming/file naming      */
void meshPhys::Refine(MAd::MeshAdapter* adapter, MAd::pMesh& mesh,
                      int array_id, int dim, std::string outMeshFile)   
{
  std::string array_name =  pntData[array_id].getName(); 
  if (dim > 1)
  {
    adapter->registerVData(array_name, fold(pntData[array_id].getData(), dim));
    // check if data is correctly attached
    std::vector<std::vector<double>> preAMRData;
    adapter->getMeshVData(array_name, &preAMRData);
    std::cout << "preAMRData[1][1] = " << preAMRData[1][1] << std::endl;
    std::cout << "Size of vector = " << preAMRData.size() << std::endl; 

    adapter->setSwapMinImproveRatio(1);

    // Output situation before refinement
    std::cout << "Statistics before refinement: " << std::endl;
    adapter->printStatistics(std::cout);
  
    // Optimize
    std::cout << "Refining the mesh ..." << std::endl;
    adapter->run();
    std::cout << "Statistics after refinement: " << std::endl;
    adapter->printStatistics(std::cout);

    // running laplacian smoothing
    std::cout << "Optimizing the mesh" << std::endl;


    for (int i = 0; i < 10; ++i)
    {
      adapter->LaplaceSmoothing();
      adapter->splitLongestEdges();
      adapter->removeSlivers();
      adapter->optimiseEdgeLength();
      adapter->optimiseElementShape();    
    }
     /* there is a fast laplace smoothing function available where instead of computing 
       the optimal position it uses the cavity center. The center is the initial position
       passed to the routine "computeOptimalLocation" before it calculates optimal. The 
       optimal way is run by default. If we modify line 714/718 in AdaptInterface.cc to
       laplOp->runFast we can use the fast method */  

 
    // Outputs final mesh
    std::cout << "Statistics after optimization: " << std::endl;
    adapter->printStatistics(std::cout);
    
    // unclassifying boundary elements for proper output
    mesh->unclassify_grid_boundaries();
    // writing refined mesh to file in msh format 
    MAd::M_writeMsh(mesh, (char*) &outMeshFile[0u], 2);
   
    // get data after refinement
    std::vector<std::vector<double>> postAMRData;
    adapter->getMeshVData(array_name, &postAMRData);
    std::cout << "postAMRData[1][1] = " << postAMRData[1][1] << std::endl;
    std::cout << "Size of vector = " << postAMRData.size() << std::endl; 

    // convert refined msh file back to vtk
    GModel* trgGModel;
    trgGModel = new GModel("refined"); 
    trgGModel->readMSH((char*) &outMeshFile[0u]);
    
    // optimizing mesh with netgen
    //std::cout << "optimize return: " << trgGModel->optimizeMesh("Netgen") << std::endl;
    //std::cout << "refine return: " << trgGModel->refineMesh(1) << std::endl;



    trgGModel->writeVTK(trim_fname(outMeshFile,".vtk"), false, true); // binary=false, saveall=true


    // write physical quantities to vtk file
    vtkAnalyzer* trgVTK;
    // TODO: Seg fault if name is like refined_smooth.vtk
    trgVTK = new vtkAnalyzer((char*) &(trim_fname(outMeshFile,".vtk"))[0u]);

    trgVTK->read();
    std::vector<double> postAMRDatas = flatten(postAMRData);
    trgVTK->setPointDataArray(&array_name[0u], dim, postAMRDatas);
    trgVTK->report();
    trgVTK->write((char*) &(trim_fname(outMeshFile, "_solution.vtu"))[0u]);
    trgVTK->writeMSH(trim_fname(outMeshFile, "_solution.msh"));
    delete trgGModel;
    delete trgVTK;
  }
  else if (dim == 1)
  {  
    adapter->registerData(array_name, pntData[array_id].getData());
    // check if data is correctly attached
    std::vector<double> preAMRData;
    adapter->getMeshData(array_name, &preAMRData);
    std::cout << "preAMRData[1]  = " << preAMRData[1] << std::endl;
    std::cout << "Size of vector = " << preAMRData.size() << std::endl; 

    // Output situation before optimization
    std::cout << "Statistics before optimization: " << std::endl;
    adapter->printStatistics(std::cout);
  
    // Optimize
    std::cout << "Optimizing the mesh ..." << std::endl;
    adapter->run();
    adapter->removeSlivers();

    // Outputs final mesh
    std::cout << "Statistics after optimization: " << std::endl;
    adapter->printStatistics(std::cout);

    // unclassifying boundary elements for proper output
    mesh->unclassify_grid_boundaries();
    
    // writing refined mesh to file in msh format 
    MAd::M_writeMsh(mesh, (char*) &outMeshFile[0u], 2);
   
    // get data after refinement
    std::vector<double> postAMRData;
    adapter->getMeshData(array_name, &postAMRData);
    std::cout << "postAMRData[1] = " << postAMRData[1] << std::endl;
    std::cout << "Size of vector = " << postAMRData.size() << std::endl; 

    // convert refined msh file back to vtk
    GModel* trgGModel;
    trgGModel = new GModel("refined"); 
    trgGModel->readMSH((char*) &outMeshFile[0u]);
    trgGModel->writeVTK(trim_fname(outMeshFile,".vtk"), false, true); // binary=false, saveall=true


    // write physical quantities to vtk file
    vtkAnalyzer* trgVTK;
    trgVTK = new vtkAnalyzer((char*) &(trim_fname(outMeshFile,".vtk"))[0u]);
    trgVTK->read();
    trgVTK->setPointDataArray(&array_name[0u], 1, postAMRData);
    trgVTK->report();
    trgVTK->write((char*) &(trim_fname(outMeshFile, "_solution.vtu"))[0u]);
    trgVTK->writeMSH(trim_fname(outMeshFile, "_solution.msh"));
    delete trgGModel;
    delete trgVTK;
  }
  else
  {
    std::cout << "Dimension of data must be >= 1!" << std::endl;
    exit(1);
  }
}


//-----------------------------------------------------------------------------------//
