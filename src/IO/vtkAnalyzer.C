#include <algorithm>

#include "vtkAnalyzer.H"

#include "baseInterp.H"

// TODO: We shouldn't be returning double arrays declared in the function
//       The stack is restored after the function's scope, so the addresses
//       may point to garbage. Either allocate them on the heap or use an stl container
//       Things work fine now, but we're just getting lucky with sequential function calls

void vtkAnalyzer::read() 
{
    extension = vtksys::SystemTools::GetFilenameLastExtension(xmlFileName);
    // Dispatch based on the file extension
    if (extension == ".vtu")
    {
      //unsGridReader = ReadAnXMLFile<vtkXMLUnstructuredGridReader> (xmlFileName);
      unsGridReader =
         vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
      unsGridReader->SetFileName(xmlFileName);
      unsGridReader->Update();
      unsGridReader->GetOutput()->Register(unsGridReader);
      dataSet = vtkUnstructuredGrid::SafeDownCast(unsGridReader->GetOutput());
    }
    else if (extension == ".vtp")
    {
      //dataSet = ReadAnXMLFile<vtkXMLPolyDataReader> (xmlFileName);
    }
    else if (extension == ".vts")
    {
      //dataSet = ReadAnXMLFile<vtkXMLStructuredGridReader> (xmlFileName);
    }
    else if (extension == ".vtr")
    {
      //dataSet = ReadAnXMLFile<vtkXMLRectilinearGridReader> (xmlFileName);
    }
   

    else if (extension == ".vti")
    {
      // adding support for image file
      imageDataSetReader = 
        vtkSmartPointer<vtkXMLImageDataReader>::New();
      imageDataSetReader->SetFileName(xmlFileName);
      imageDataSetReader->Update();
      imageDataSetReader->GetOutput()->Register(imageDataSetReader);
      dataSet = vtkDataSet::SafeDownCast(imageDataSetReader->GetOutput());
    }
    // NOTE: There seems to be an issue with reading array data in .vtk format
    else if (extension == ".vtk")
    {
      //dataSetReader = ReadAnXMLFile<vtkDataSetReader> (xmlFileName);
      dataSetReader =
         vtkSmartPointer<vtkDataSetReader>::New();
      dataSetReader->SetFileName(xmlFileName);
      dataSetReader->Update();
      dataSetReader->GetOutput()->Register(dataSetReader);
      dataSet = vtkDataSet::SafeDownCast(dataSetReader->GetOutput());
    }
    else {
      std::cerr << " Unknown file format for " << xmlFileName << std::endl;
    }
}

void vtkAnalyzer::write(std::string outXMLFileName)
{

    // Dispatch based on the file extension
    if (extension == ".vtu")
    {
      if (unsGridReader) {
        WriteAnXMLFile<vtkXMLUnstructuredGridWriter> (outXMLFileName.c_str(), unsGridReader->GetOutput());
      }
    }
    else if (extension == ".vtp")
    {
      //dataSet = WriteAnXMLFile<vtkXMLPolyDataReader> (outXMLFileName);
    }
    else if (extension == ".vts")
    {
      //dataSet = ReadAnXMLFile<vtkXMLStructuredGridReader> (xmlFileName);
    }
    else if (extension == ".vtr")
    {
      //dataSet = ReadAnXMLFile<vtkXMLRectilinearGridReader> (xmlFileName);
    }
    else if (extension == ".vti")
    {
      //dataSet = ReadAnXMLFile<vtkXMLImageDataReader> (xmlFileName);
    }
    else if (extension == ".vtk")
    {
      if (dataSetReader) {
        WriteAnXMLFile<vtkXMLUnstructuredGridWriter> (outXMLFileName.c_str(), dataSetReader->GetOutput());
      }
    }
    else {
      std::cerr << " Unknown file format for " << xmlFileName << std::endl;
    }
}
int vtkAnalyzer::getNumberOfPoints() 
{
   numberOfPoints = dataSet->GetNumberOfPoints();
   return numberOfPoints;
}

int vtkAnalyzer::getNumberOfCells() 
{
   numberOfCells = dataSet->GetNumberOfCells();
   return numberOfCells;
}

// check for named array in vtk 
int vtkAnalyzer::IsArrayName(std::string name)
{
  vtkPointData* pd = dataSet->GetPointData();
  if (pd->GetNumberOfArrays()) {
    for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
      std::string curr_name = (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL");
      if (!name.compare(curr_name)) {
        return i;
      }
    }
    // fall through to exit
    std::cout << "Invalid Species Name" << std::endl;
    std::cout << "Valid options are:" << std::endl;
    for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
      std::cout << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")
                << std::endl;
    }
    exit(1);
  }
  return -1;
}

void vtkAnalyzer::report() 
{
   // populate if not yet
   if (!dataSet)
      read();

   // Generate a report
   std::cout << "Processing the file ..... " << std::endl;
   std::cout << xmlFileName << std::endl
      << " contains a " 
      << dataSet->GetClassName()
      <<  " that has " << getNumberOfCells() << " cells"
      << " and " << getNumberOfPoints() << " points." << std::endl;

   CellContainer cellMap;
   for (int i = 0; i < numberOfCells; i++)
   {
   cellMap[dataSet->GetCellType(i)]++;
   }

   CellContainer::const_iterator it = cellMap.begin();
   while (it != cellMap.end())
   {
   std::cout << "\tCell type "
  << vtkCellTypes::GetClassNameFromTypeId(it->first)
  << " occurs " << it->second << " times." << std::endl;
   ++it;
   }

   // Now check for point data
   vtkPointData *pd = dataSet->GetPointData();
   if (pd)
   {
     std::cout << " contains point data with "
          << pd->GetNumberOfArrays()
          << " arrays." << std::endl;
     for (int i = 0; i < pd->GetNumberOfArrays(); i++)
     {
       std::cout << "\tArray " << i
            << " is named "
            << (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL") ;
       vtkDataArray* da = pd->GetArray(i);
       std::cout << " with " << da->GetNumberOfTuples() 
                 << " values. " << std::endl;
     }
   }

   // Now check for cell data
   vtkCellData *cd = dataSet->GetCellData();
   if (cd)
   {
     std::cout << " contains cell data with "
          << cd->GetNumberOfArrays()
          << " arrays." << std::endl;
     for (int i = 0; i < cd->GetNumberOfArrays(); i++)
      {
       std::cout << "\tArray " << i
            << " is named "
            << (cd->GetArrayName(i) ? cd->GetArrayName(i) : "NULL") ;
            vtkDataArray* da = cd->GetArray(i);
            std::cout << " with " << da->GetNumberOfTuples() 
                      << " values. " << std::endl;
      }
   }
   // Now check for field data
   if (dataSet->GetFieldData())
   {
      std::cout << " contains field data with "
           << dataSet->GetFieldData()->GetNumberOfArrays()
           << " arrays." << std::endl;
      for (int i = 0; i < dataSet->GetFieldData()->GetNumberOfArrays(); i++)
      {
        std::cout << "\tArray " << i
             << " is named " << dataSet->GetFieldData()->GetArray(i)->GetName();
            vtkDataArray* da = dataSet->GetFieldData()->GetArray(i);
            std::cout << " with " << da->GetNumberOfTuples() 
                      << " values. " << std::endl;
      }
   }
}

int vtkAnalyzer::getNumberOfPointData()
{
   pointData = dataSet->GetPointData();
   if (pointData)
   {
    numberOfPointData = pointData->GetNumberOfArrays();
   }
   return numberOfPointData;
}

int vtkAnalyzer::getNumberOfCellData()
{
   cellData = dataSet->GetCellData();
   if (cellData)
   {
    numberOfCellData = cellData->GetNumberOfArrays();
   }
   return numberOfCellData;
}

// if cell face belongs to only 1 cell, it is a surface element
std::multimap<int, std::vector<int> > vtkAnalyzer::findBoundaryFaces()
{
  int numCells = getNumberOfCells();
  vtkCell* cell;
  vtkCell* face;

  int npts;
  std::multimap<int, std::vector<int> > boundaries;

  for (int i = 0; i < numCells; ++i)
  {

    cell = dataSet->GetCell(i);
    if (cell->GetCellDimension() == 3)
    {

      int numFaces = cell->GetNumberOfFaces();
      for (int j = 0; j < numFaces; ++j)
      {
        face = cell->GetFace(j);
        vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> facePntIds = vtkSmartPointer<vtkIdList>::New();
        facePntIds = face->GetPointIds(); 
        dataSet->GetCellNeighbors(i,facePntIds,cellIds.GetPointer());
        if (cellIds->GetNumberOfIds() <= 0 )
        {

          npts = face->GetNumberOfPoints();
          std::vector<int> ptIds(npts);
          for (int k = 0; k < npts; ++k)
          {
              ptIds[k] = face->GetPointId(k); 
          }
          boundaries.insert(std::pair<int,std::vector<int> > (i,ptIds)); 
        }
      }
    }
    else if (cell->GetCellDimension() == 2)
    {
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      vtkSmartPointer<vtkIdList> facePntIds = vtkSmartPointer<vtkIdList>::New();
      facePntIds = cell->GetPointIds(); 
      dataSet->GetCellNeighbors(i,facePntIds,cellIds.GetPointer());
      if (cellIds->GetNumberOfIds() <= 1 )
      {
        npts = facePntIds->GetNumberOfIds();
        std::vector<int> ptIds(npts);
        for (int k = 0; k < npts; ++k)
        {
            ptIds[k] = cell->GetPointId(k); 
        }
        boundaries.insert(std::pair<int,std::vector<int> > (i,ptIds)); 
      }
      
    }
  }

  return boundaries; 
}


 
std::vector<std::vector<double*>> vtkAnalyzer::getSurfaceTriElements(int& numComponent)
{
  getNumberOfCells();
  std::vector<std::vector<double*>> coords;
  for (int j = 0; j < numberOfCells; ++j)
  {
    if (dataSet->GetCellType(j) == VTK_TRIANGLE)
    { 
      coords.push_back(getCellCoords(j, numComponent));
    } 
  }
  return coords;
}

void vtkAnalyzer::writeSurfaceTriElements(std::string fname)
{
    
  std::ofstream vtk(fname);

  if (!vtk.good())
  {
    std::cout << "Error opening: " << fname << std::endl;
    exit(1);
  }

  int numComponent;
  std::vector<std::vector<double*>> surfTri = getSurfaceTriElements(numComponent);
  vtk << "# vtk DataFile Version 2.0" << std::endl 
      << "tmp surf mesh" << std::endl
      << "ASCII" << std::endl
      << "DATASET UNSTRUCTURED_GRID" << std::endl 
      << "POINTS " << surfTri.size()*3 << " double" << std::endl;

  for (int i = 0; i < surfTri.size(); ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        vtk << surfTri[i][j][k] << " ";
      }
      delete surfTri[i][j];
      vtk << std::endl; 
    }
  }
  vtk << "CELLS " << surfTri.size() << " " << surfTri.size()*4 << std::endl;
  for (int i = 0; i < surfTri.size(); ++i)
  {
    vtk << 3 << " ";
    for (int j = 0; j < 3; ++j)
      vtk << i*3 +j << " ";
    vtk << std::endl;
  }
  vtk << "CELL_TYPES " << surfTri.size() << std::endl;
  for (int i = 0; i < surfTri.size(); ++i)
    vtk << 5 << std::endl; 
 
}


double* vtkAnalyzer::getPointCoords(int pntId)
{
   double* pntCoords;
   if (pntId < dataSet->GetNumberOfPoints()) {
     pntCoords = dataSet->GetPoint(pntId);
   } else {
     std::cerr << "Point ID is out of range!" << std::endl;
     exit(1);
   }
   return pntCoords;
}

// returns coordinates of member points in cell with ID
// call delete when finished with return
std::vector<double* > vtkAnalyzer::getCellCoords(int cellId, int& numComponent)
{
  std::vector<double *> cellCoords;
  if (cellId < dataSet->GetNumberOfCells()) {
    vtkIdList* point_ids = dataSet->GetCell(cellId)->GetPointIds();
    numComponent = point_ids->GetNumberOfIds();
    cellCoords.resize(numComponent);// = new double*[numComponent];
    for (int i = 0; i < numComponent; ++i) {
      cellCoords[i] = new double[3];
      dataSet->GetPoint(point_ids->GetId(i), cellCoords[i]);  
    }
  }
  else {
    std::cerr << "Cell ID is out of range!" << std::endl;
    exit(2);
  }
  return cellCoords;
}

// get number of non triangular/line/point elements
int vtkAnalyzer::getNumberOfNonTri()
{
  getNumberOfCells(); 
  int num = 0;
  for (int j = 0; j < numberOfCells; ++j)
  {
    if (dataSet->GetCellType(j) < VTK_TRIANGLE)
      num++;
  }
  return num;
}

// returns coordinates of all points in mesh
std::vector<double> vtkAnalyzer::getAllPointCoords(int nDim)
{
  int num_points = getNumberOfPoints();
  std::vector<double> VolPointCoords(num_points*nDim);
  for (int i = 0; i < num_points; ++i) {
    double* pntcrds = getPointCoords(i);
    VolPointCoords[i*nDim] = pntcrds[0];
    VolPointCoords[i*nDim+1] = pntcrds[1];
    VolPointCoords[i*nDim+2] = pntcrds[2];
    
  }
return VolPointCoords;
}

// find min extent
double vtkAnalyzer::getMinExtent(int nDim, const std::vector<double>& pntCrds)
{
  // creating vector for x, y and z
  std::vector<double> xCrds, yCrds,zCrds;
  for (int i = 0; i < pntCrds.size(); i+=nDim) {
    xCrds.push_back(pntCrds[i]);
    yCrds.push_back(pntCrds[i+1]);
    zCrds.push_back(pntCrds[i+2]);
  }
  
  // using minmax_element from stl
  auto Xminmax = std::minmax_element(xCrds.begin(), xCrds.end());
  auto Yminmax = std::minmax_element(yCrds.begin(), yCrds.end());
  auto Zminmax = std::minmax_element(zCrds.begin(), zCrds.end());

  // populating extent vector
  std::vector<double> extents;
  extents.push_back(*Xminmax.second - *Xminmax.first);
  extents.push_back(*Yminmax.second - *Yminmax.first);
  extents.push_back(*Zminmax.second - *Zminmax.first);
  // returning value of min element from extent vector
  return *std::min_element(extents.begin(), extents.end()); 
} 


// returns centers of all cells
std::vector<double>
vtkAnalyzer::getCellCenters(int& numComponent)
{

  
  int num_cells = getNumberOfCells();
  // enforce that cell types are not lines etc
  std::vector<double> cellCenters;
  for (int i = 0; i < num_cells; ++i) {
   // if (dataSet->GetCellType(i) != VTK_TRIANGLE)
   //   deleteCell(i);
   // else {
    
    std::vector<double *> cellCoords = getCellCoords(i, numComponent);
    double x, y,z;
    x = y = z = 0;
    for (int j = 0; j < numComponent ; ++j) {
      x += cellCoords[j][0]/numComponent;
      y += cellCoords[j][1]/numComponent;
      z += cellCoords[j][2]/numComponent;
      delete cellCoords[j];
    }
    cellCenters.push_back(x);
    cellCenters.push_back(y);
    cellCenters.push_back(z);
   // }
  }
  return cellCenters;
}   

// search for pnt id in each cell and return common
std::vector<int> vtkAnalyzer::getCellsWithPoint(int pnt)
{
  getNumberOfCells();
    
  std::vector<int> commonCells;
  for (int i = 0; i < numberOfCells; ++i)
  {
    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New(); 
    //point_ids = dataSet->GetCell(i)->GetPointIds();
    dataSet->GetCellPoints(i, point_ids);
    int numComponent = point_ids->GetNumberOfIds();
    for (int j = 0; j < numComponent; ++j)
    {
      if (point_ids->GetId(j) == pnt)
      {
        commonCells.push_back(i);
        std::cout << "point is in cell: " << i  << std::endl;
      }
    }
  }
  return commonCells;  
}
    

// interpolate point data from 3D mesh in neighborhoods of
// cell centers of planar mesh to those centers.
std::vector<std::vector<double>>
vtkAnalyzer::getInterpData(int nDim, int num_neighbors, int numComponent, int numTuple,
                           std::vector<std::vector<double>>& volDataMat,
                           std::vector<double>& PlaneCellCenters,
                           std::vector<double>& VolPointCoords, double tol)
{
  //std::vector<double> VolPointCoords = getAllPointCoords(nDim);
  int num_vol_points = getNumberOfPoints(); 
  int num_interp_points = PlaneCellCenters.size()/nDim; 

  basicInterpolant* VolPointInterp = 
    new basicInterpolant(nDim, num_vol_points, num_neighbors, VolPointCoords);

  std::vector<std::vector<double>> interpData(numComponent);  
  for (int j = 0; j < numComponent; ++j) {  
    std::vector<double > volData(numTuple);
    for (int i = 0; i < numTuple; ++i) {
      volData[i] = volDataMat[i][j];
    }   
    VolPointInterp->interpolate(num_interp_points, 
                               PlaneCellCenters, volData, interpData[j], tol, 0);
  }
  delete VolPointInterp;
  return interpData;
}

// interpolate point data from 3D mesh in neighborhoods of
// cell centers of planar mesh to those centers for cases with spheres
// numTuple should = VolPointCoords.size()/ndim
std::vector<std::vector<double>>
vtkAnalyzer::getInterpData(int nDim, int num_neighbors, int numComponent, int numTuple,
                           std::vector<std::vector<double>>& volDataMat,
                           std::vector<double>& PlaneCellCenters,
                           std::vector<double>& VolPointCoords,
                           std::vector<sphere>& spheres, std::vector<double>& maskData, double tol)
{
  //std::vector<double> VolPointCoords = getAllPointCoords(nDim);
  int num_vol_points = getNumberOfPoints(); 
  int num_interp_points = PlaneCellCenters.size()/nDim; 

  basicInterpolant* VolPointInterp = 
    new basicInterpolant(nDim, num_vol_points, num_neighbors, VolPointCoords);

  std::vector<std::vector<double>> interpData(numComponent);  
  for (int j = 0; j < numComponent; ++j) {   
    std::vector<double > volData(numTuple);
    for (int i = 0; i < numTuple; ++i) {
      // volData already 0 in sphere from RocLB - bool in_sphere = false;
      volData[i] = volDataMat[i][j];
      
    }   
    VolPointInterp->interpolate(num_interp_points,         
                                PlaneCellCenters, spheres, maskData, volData, interpData[j], tol, 0);
  }
  delete VolPointInterp;
  return interpData;
}

// no spheres, write_coords=0
void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default,
                                  double poisson_dom_default, double T, double R,
                                  std::ostream& outputStream, int numNonTri,
                                  const double conc_conv)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  }

  outputStream << std::left << std::setw(10) << "id" 
               << std::left << std::setw(16) << "rho"
               << std::left << std::setw(16) << "E"
               << std::left << std::setw(16) << "V" << std::endl << std::endl;
  for (int i = 0; i < interpData[0].size(); ++i) {
    if (i > numNonTri - 1) {
      outputStream << std::left << std::setw(10) << i-numNonTri << std::left << std::setw(16); 
      for (int j = 0; j < interpData.size(); ++j) {
        outputStream << interpData[j][i] << std::left << std::setw(16);
      }
      double G = interpData[0][i]*R*T*(1 - Mc/M) * conc_conv; //last term converts rho to SI units
      // if V not given, get from G and E
      if (poisson_dom_default == -1) {
        double V = youngs_dom_default/(2*G) - 1;
      // correcting for negative 0 output
        outputStream << std::left << std::setw(16) << youngs_dom_default 
                     << std::left << std::setw(16) << (V == 0.0 ? fabs(V): V) << std::endl;
      } 
      // if E not given, get from G and V
      else if (youngs_dom_default == -1) {
        double E = 2.0*G*(1.0+poisson_dom_default);
        outputStream << std::left << std::setw(16) << E 
                     << std::left << std::setw(16) << poisson_dom_default << std::endl;
      }
      else {
        outputStream << std::left << std::setw(16) << youngs_dom_default
                     << std::left << std::setw(16) << poisson_dom_default << std::endl;
      } 
    }
  }   
}

// no spheres and coord switch
void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default,
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::ostream& outputStream, bool writeCoord, int numNonTri,
                                  const double conc_conv)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  }
  if (!writeCoord)
    writeInterpData(interpData, Mc, M, youngs_dom_default,
                    poisson_dom_default, T, R,
                    outputStream, numNonTri, conc_conv);
  else {
    outputStream << std::left << std::setw(10) << "id" 
                 << std::left << std::setw(16) << "X"
                 << std::left << std::setw(16) << "Y"
                 << std::left << std::setw(16) << "Z" 
                 << std::left << std::setw(16) << "rho"   
                 << std::left << std::setw(16) << "E"
                 << std::left << std::setw(16) << "V" << std::endl << std::endl;
    //double R = .000008314;
    //double T = 300.0; 
    for (int i = 0; i < interpData[0].size(); ++i) {
      if (i > numNonTri-1) {
        outputStream << std::left << std::setw(10) << i-numNonTri
                     << std::left << std::setw(16) << PlaneCellCenters[i*nDim] 
                     << std::left << std::setw(16) << PlaneCellCenters[i*nDim+1] 
                     << std::left << std::setw(16) << PlaneCellCenters[i*nDim+2] ;
        for (int j = 0; j < interpData.size(); ++j) {
          outputStream << std::left << std::setw(16) << interpData[j][i] * conc_conv
                       << std::left << std::setw(16);
        }  
        double G = interpData[0][i]*R*T*(1 - Mc/M) * conc_conv;  //last term converts unit to SI
        // if V not given, get from G and E
        if (poisson_dom_default == -1) {
          double V = youngs_dom_default/(2*G) - 1;
        // correcting for negative 0 output
          outputStream << std::left << std::setw(16) << youngs_dom_default 
                       << std::left << std::setw(16) << (V == 0.0 ? fabs(V): V) << std::endl;
        } 
        // if E not given, get from G and V
        else if (youngs_dom_default == -1) {
          double E = 2.0*G*(1.0+poisson_dom_default);
          outputStream << std::left << std::setw(16) << E 
                       << std::left << std::setw(16) << poisson_dom_default << std::endl;
        }
        else {
          outputStream << std::left << std::setw(16) << youngs_dom_default
                       << std::left << std::setw(16) << poisson_dom_default << std::endl;
        } 
      }
    }
  }    
}

// has spheres, write_coords=0
void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default,
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::vector<sphere>& spheres,
                                  std::vector<std::string>& mat_sphere_names,
                                  std::vector<std::string>& material_names,
                                  std::vector<double>& youngs_inc_default,
                                  std::vector<double>& shear_inc_default,
                                  std::vector<double>& poisson_inc_default,
                                  std::ostream& outputStream, int numNonTri,
                                  const double conc_conv)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  }
  
  outputStream << std::left << std::setw(10) << "id" 
               << std::left << std::setw(16) << "rho"
               << std::left << std::setw(16) << "E"
               << std::left << std::setw(16) << "V" << std::endl << std::endl;
  
  //double R = .000008314;
  //double T = 300.0; 
  for (int i = 0; i < interpData[0].size(); ++i) {
    if (i > numNonTri - 1) {
      outputStream << std::left << std::setw(10) << i-numNonTri << std::left << std::setw(16); 
      bool in_sphere = false;
      int k;
      for (int j = 0; j < interpData.size(); ++j) {
        // checking if plane points are in sphere
        std::vector<double> point;
        point.push_back(PlaneCellCenters[i*nDim]);
        point.push_back(PlaneCellCenters[i*nDim+1]);
        point.push_back(PlaneCellCenters[i*nDim+2]);
        for (k = 0; k < spheres.size(); ++k) {
          in_sphere = spheres[k].in_sphere(point);
          if(in_sphere) {
            break;
          }
        }
      //if (in_sphere)
      //  outputStream << 0.0 << std::left << std::setw(16);
      //else      
        outputStream << interpData[j][i] * conc_conv << std::left << std::setw(16);
      }
      if (in_sphere) {
        for (int i = 0; i < material_names.size(); ++i) {
          if (!material_names[i].compare(mat_sphere_names[k])) {
            // if V not given, get from G and E
            if (poisson_inc_default[i] == -1) {
              outputStream << std::left << std::setw(16) << youngs_inc_default[i]
                           << std::left << std::setw(16)
                           << (youngs_inc_default[i])/(2*shear_inc_default[i]) - 1
                           << std::endl;
            }
            // if E not given, get from G and V
            else if (youngs_inc_default[i] == -1) {
              outputStream << std::left << std::setw(16)
                           << 2.0*shear_inc_default[i]*(1.0 + poisson_inc_default[i])
                           << std::left << std::setw(16) << poisson_inc_default[i] << std::endl;
            }
            // if E and V given, print
            else {
              outputStream << std::left << std::setw(16)
                           << youngs_inc_default[i]
                           << std::left << std::setw(16)
                           << poisson_inc_default[i] << std::endl;
            }
          }
        }
      }
      else {
        double G = interpData[0][i]*R*T*(1 - Mc/M) * conc_conv; //last term converts unit to SI
        // if V not given, get from G and E
        if (poisson_dom_default == -1) {
        double V = youngs_dom_default/(2*G) - 1;
        // correcting for negative 0 output
        outputStream << std::left << std::setw(16) << youngs_dom_default
                     << std::left << std::setw(16) << (V == 0.0 ? fabs(V): V) << std::endl;
        }
        // if E not given, get from G and V
        else if (youngs_dom_default == -1) {
        double E = 2.0*G*(1.0+poisson_dom_default);
        outputStream << std::left << std::setw(16) << E
                     << std::left << std::setw(16) << poisson_dom_default << std::endl;
        }
        else {
          outputStream << std::left << std::setw(16) << youngs_dom_default
                       << std::left << std::setw(16) << poisson_dom_default << std::endl;
        }
      }
    }
  }   
}

// has spheres and coord switch
void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default, 
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::vector<sphere>& spheres,
                                  std::vector<std::string>& mat_sphere_names,
                                  std::vector<std::string>& material_names,
                                  std::vector<double>& youngs_inc_default,
                                  std::vector<double>& shear_inc_default,
                                  std::vector<double>& poisson_inc_default,
                                  std::ostream& outputStream, bool writeCoord, int numNonTri,
                                  const double conc_conv)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  } 
  if (!writeCoord)
    writeInterpData(interpData, Mc, M, youngs_dom_default, 
                    poisson_dom_default, T, R,
                    PlaneCellCenters, nDim, spheres, 
                    mat_sphere_names, material_names, 
                    youngs_inc_default, shear_inc_default, 
                    poisson_inc_default, outputStream, numNonTri, conc_conv);
  else {
    outputStream << std::left << std::setw(10) << "id" 
                 << std::left << std::setw(16) << "X"
                 << std::left << std::setw(16) << "Y"
                 << std::left << std::setw(16) << "Z" 
                 << std::left << std::setw(16) << "rho"   
                 << std::left << std::setw(16) << "E"
                 << std::left << std::setw(16) << "V" << std::endl << std::endl;
    
    for (int i = 0; i < interpData[0].size(); ++i) {
      if (i > numNonTri-1) {
        outputStream << std::left << std::setw(10) << i - numNonTri  
                     << std::left << std::setw(16) << PlaneCellCenters[i*nDim] 
                     << std::left << std::setw(16) << PlaneCellCenters[i*nDim+1] 
                     << std::left << std::setw(16) << PlaneCellCenters[i*nDim+2] ;
                          //<< std::left << std::setw(16);
        bool in_sphere = false;
        int k; // sphere-material identifier
        for (int j = 0; j < interpData.size(); ++j) {
          // checking if plane points are in sphere
          std::vector<double> point;
          point.push_back(PlaneCellCenters[i*nDim]);
          point.push_back(PlaneCellCenters[i*nDim+1]);
          point.push_back(PlaneCellCenters[i*nDim+2]);
          for (k = 0; k < spheres.size(); ++k) {
            in_sphere = spheres[k].in_sphere(point);
            if(in_sphere) {
              break;
            }
          }
        // if point in plane is in sphere, crosslink = 0
        //if (in_sphere) 
        //  outputStream << std::left << std::setw(16) << 0.0 
        //               << std::left << std::setw(16);
        //else      
          outputStream << std::left << std::setw(16) << interpData[j][i] * conc_conv 
                       << std::left << std::setw(16);
        } 
        if (in_sphere) {
          for (int i = 0; i < material_names.size(); ++i) { 
            if (!material_names[i].compare(mat_sphere_names[k])) {
              // if V not given, get from G and E
              if (poisson_inc_default[i] == -1) {
                outputStream << std::left << std::setw(16) << youngs_inc_default[i] 
                             << std::left << std::setw(16) 
                             << (youngs_inc_default[i])/(2*shear_inc_default[i]) -1
                             << std::endl;
              }
              // if E not given, get from G and V
              else if (youngs_inc_default[i] == -1) {
                outputStream << std::left << std::setw(16) 
                             << 2.0*shear_inc_default[i]*(1.0 + poisson_inc_default[i])
                             << std::left << std::setw(16) << poisson_inc_default[i] << std::endl;
              }   
              // if E and V given, print
              else {
                outputStream << std::left << std::setw(16) 
                             << youngs_inc_default[i]
                             << std::left << std::setw(16) 
                             << poisson_inc_default[i] << std::endl;
              }
            }
          }
        }
        else {
          double G = interpData[0][i]*R*T*(1 - Mc/M) * conc_conv; //last term converts unit to SI
          // if V not given, get from G and E
          if (poisson_dom_default == -1) {
            double V = youngs_dom_default/(2*G) - 1;
          // correcting for negative 0 output
            outputStream << std::left << std::setw(16) << youngs_dom_default 
                         << std::left << std::setw(16) << (V == 0.0 ? fabs(V): V) << std::endl;
          } 
          // if E not given, get from G and V
          else if (youngs_dom_default == -1) {
            double E = 2.0*G*(1.0+poisson_dom_default);
            outputStream << std::left << std::setw(16) << E 
                        << std::left << std::setw(16) << poisson_dom_default << std::endl;
          }
          else {
            outputStream << std::left << std::setw(16) << youngs_dom_default
                         << std::left << std::setw(16) << poisson_dom_default << std::endl;
          } 
        }
      }
    }
  }
}



void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default, 
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::string filename, bool writeCoord, int numNonTri,
                                  const double conc_conv)
{
  std::ofstream outputStream(filename.c_str());
  if(!outputStream.good()) {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
  writeInterpData(interpData, Mc, M, youngs_dom_default, 
                  poisson_dom_default, T, R, PlaneCellCenters, nDim, 
                  outputStream, writeCoord, numNonTri, conc_conv);
}                                 


void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                       double Mc, double M, double youngs_dom_default, 
                       double poisson_dom_default, double T, double R,
                       const std::vector<double>& PlaneCellCenters, int nDim,
                       std::vector<sphere>& spheres,
                       std::vector<std::string>& mat_sphere_names,
                       std::vector<std::string>& material_names,
                       std::vector<double>& youngs_inc_default,
                       std::vector<double>& shear_inc_default,
                       std::vector<double>& poisson_inc_default,
                       std::string filename, bool writeCoord, int numNonTri,
                       const double conc_conv)
{
  std::ofstream outputStream(filename.c_str());
  if(!outputStream.good()) {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
  writeInterpData(interpData, Mc, M, youngs_dom_default, poisson_dom_default, T, R,
                  PlaneCellCenters, nDim, spheres, mat_sphere_names, material_names, 
                  youngs_inc_default, shear_inc_default, poisson_inc_default, 
                  outputStream, writeCoord, numNonTri, conc_conv);
}                                 

void vtkAnalyzer::writeCSV(char* fname, std::vector<double> slnVec)
{
   std::ofstream csvFile;
   csvFile.open(fname);
   csvFile << "x, y, z, T" << std::endl;
   for (int iPnt=0; iPnt<dataSet->GetNumberOfPoints(); iPnt++){
   double* pntCoords = getPointCoords(iPnt);
   csvFile << pntCoords[0] << ", " 
     << pntCoords[1] << ", "
     << pntCoords[2] << ", "
     << slnVec[iPnt]
     << std::endl;          
   }
   csvFile.close();
} 

int vtkAnalyzer::getPointDataArray(int id, std::vector<std::vector<double> > &pntData, 
                                   int &numTuple, int &numComponent)
{
   // TODO: This function is slow
   // check data exists
  // vtkDoubleArray* data;
  if (!pointData)
    pointData = dataSet->GetPointData();
   // populate user's array
  if (pointData) {
    vtkDataArray* da = pointData->GetArray(id);
    if (da) {
      numComponent = da->GetNumberOfComponents();
      numTuple = da->GetNumberOfTuples();
      pntData.resize(numTuple);
      for (int iTup=0; iTup<numTuple; iTup++){
        double* comps = new double[numComponent];
        // EDIT: causing indexing issues
        //       shouldn't resize if pushing back //  pntData[iTup].resize(numComponent);
        da->GetTuple(iTup, comps);
        for (int iComp=0; iComp<numComponent; iComp++){
          pntData[iTup].push_back(comps[iComp]);
          }
      }
        return(numComponent*numTuple);
      } 
      else {
        return(0);
      }
  } 
  else {
    std::cerr << "There is no point data exisiting in" 
              << xmlFileName << std::endl;
    return(0);
  }
}

int vtkAnalyzer::getCellDataArray(int id, std::vector<std::vector<double> > &cllData, 
                                   int &numTuple, int &numComponent)
{
   // TODO: This function is slow
   // check data exists
   // vtkDoubleArray* data;
   if (!cellData)
    cellData = dataSet->GetCellData();
   // populate user's array
   if (cellData)
   {
   vtkDataArray* da = cellData->GetArray(id);
   if (da) {
      numComponent = da->GetNumberOfComponents();
      numTuple = da->GetNumberOfTuples();
      cllData.resize(numTuple);
      for (int iTup=0; iTup<numTuple; iTup++){
    double* comps;
    cllData[iTup].resize(numComponent);
    comps = da->GetTuple(iTup);
    for (int iComp=0; iComp<numComponent; iComp++){
      cllData[iTup].push_back(comps[iComp]);
    }
      }
      return(numComponent*numTuple);
   }
   else {
      return(0);
   }
   } 
   else {
   std::cerr << "There is no point data exisiting in" 
       << xmlFileName << std::endl;
   return(0);
   }
}

// fixed error with implementation
void vtkAnalyzer::setPointDataArray(const char* name, int numComponent, 
                                    std::vector<double> &pntArrData) {
   if (!pointData)
    pointData = dataSet->GetPointData();
   if (pointData)
   { 
    vtkSmartPointer<vtkDoubleArray> da = 
      vtkSmartPointer<vtkDoubleArray>::New();
    da->SetName(name);
    da->SetNumberOfComponents(numComponent);
    for(int i=0; i<getNumberOfPoints(); i++)
    {
      double* pntData = new double[numComponent];
      for (int j = 0; j < numComponent; ++j)
        pntData[j] = pntArrData[i*numComponent + j];
      da->InsertNextTuple(pntData);
	  delete [] pntData;
    }
    pointData->SetActiveScalars(name);
    pointData->SetScalars(da);
   
  }
}

void vtkAnalyzer::setCellDataArray(const char* name, int numComponent,
                                    std::vector<double> &cellArrData)
{
  if (!cellData)
   cellData = dataSet->GetCellData();
  if (cellData)
  { 
    vtkSmartPointer<vtkDoubleArray> da = 
    vtkSmartPointer<vtkDoubleArray>::New();
    da->SetName(name);
    da->SetNumberOfComponents(numComponent);
    for(int i=0; i<getNumberOfCells(); i++)
    {
      double* cllData = new double[numComponent];
      for (int j = 0; j < numComponent; ++j)
        cllData[j] = cellArrData[i*numComponent +j];
      da->InsertNextTuple(cllData);
      cellData->SetActiveScalars(name);
      cellData->SetScalars(da);
	  delete[] cllData;
    } 
  }
}

void vtkAnalyzer::writeMSH(std::string filename)
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
  for (int i = 0; i < numberOfCells; i++)
  {
    int type_id = dataSet->GetCellType(i);
    if (!(type_id == 3 || type_id == 5 || type_id == 10))
    {
      std::cout << "Error: Only tetrahedral and triangular" 
                << " meshes can be written to gmsh format" << std::endl;
      exit(3);
    }
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
  outputStream << "$Elements" << std::endl << numberOfCells << std::endl;
  for (int i = 0; i < numberOfCells; ++i)
  {

    vtkIdList* point_ids = dataSet->GetCell(i)->GetPointIds();
    int numComponent = point_ids->GetNumberOfIds();
    outputStream << i + 1 << " ";
    switch(numComponent)
    {
      case 2:
      {
        outputStream << 1 << " " << 2 << " " << 1 << " " << 1 << " ";
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
        std::cerr << "Components in cell should be <= 4"<< std::endl;
        exit(1);
      }
    }
    for (int j = 0; j < numComponent; ++j)
       outputStream << point_ids->GetId(j) + 1 << " ";
    outputStream << std::endl;
  }
  outputStream << "$EndElements" << std::endl;
    
  // ---------------------------- write point data ------------------------- //

  if (!pointData)
    pointData = dataSet->GetPointData(); 
  if (pointData)
  {
    int num_arrays = pointData->GetNumberOfArrays();
    for (int i = 0; i < num_arrays; ++i)
    {
      vtkDataArray* da = pointData->GetArray(i);
      if(da)
      {
        int numComponent = da->GetNumberOfComponents();
        int numTuple = da->GetNumberOfTuples();
        std::string tmpname = "PointArray";
        tmpname += std::to_string(i);
        outputStream << "$NodeData" << std::endl
                     << 1 << std::endl // 1 string tag
                     << "\"" << (pointData->GetArrayName(i) ? 
                                 pointData->GetArrayName(i) : tmpname) // name of view
                     << "\"" << std::endl 
                     << 0 << std::endl // 0 real tag
                     << 3 << std::endl // 3 int tags (dt index, dim of field, number of fields)
                     << 0 << std::endl // dt index
                     << numComponent << std::endl // dim of field
                     << numTuple << std::endl; // number of fields
        for (int j = 0; j < numTuple; ++j)
        {
          double* data = da->GetTuple(j);
          outputStream << j + 1 << " ";
          for (int k = 0; k < numComponent; ++k)
          {
            outputStream << data[k] << " ";
          }
          outputStream << std::endl;
        }
        outputStream << "$EndNodeData" << std::endl;    
      }
    }
  }
  
  // -------------------------- write cell data ---------------------------- // 

  if (!cellData)
    cellData = dataSet->GetCellData();
  if (cellData)
  {
    int num_arrays = cellData->GetNumberOfArrays();
    for (int i = 0; i < num_arrays; ++i)
    {
      vtkDataArray* da = cellData->GetArray(i);
      if (da)
      {
        int numComponent = da->GetNumberOfComponents();
        int numTuple = da->GetNumberOfTuples();
        std::string tmpname = "CellArray";
        tmpname += std::to_string(i);
        outputStream << "$ElementData" << std::endl
                     << 1 << std::endl // 1 string tag
                     << "\"" << (cellData->GetArrayName(i) ? 
                                 cellData->GetArrayName(i) : tmpname) // name of view
                     << "\"" << std::endl 
                     << 0 << std::endl // 0 real tag
                     << 3 << std::endl // 3 int tags (dt index, dim of field, number of fields)
                     << 0 << std::endl // dt index
                     << numComponent << std::endl // dim of field
                     << numTuple << std::endl; // number of fields
        for (int j = 0; j < numTuple; ++j)
        {
          double* data = da->GetTuple(j);
          outputStream << j + 1 << " ";
          for (int k = 0; k < numComponent; ++k)
          {
            outputStream << data[k] << " ";
          }
          outputStream << std::endl;
        }
        outputStream << "$EndElementData" << std::endl;    
      }
    }
  }
}


void vtkAnalyzer::writeBackgroundMSH(std::string filename, const double size)
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
      outputStream << size << " ";
      outputStream << std::endl;
    //  i+=1;
    }
  }
  outputStream << "$EndElementData" << std::endl;    

}
