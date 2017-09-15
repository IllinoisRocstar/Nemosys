#include <vtkAnalyzer.H>

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
      dataSet = vtkDataSet::SafeDownCast(unsGridReader->GetOutput());
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
    else if (extension == ".vto")
    {
      //dataSet = ReadAnXMLFile<vtkXMLHyperOctreeReader> (xmlFileName);
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

void vtkAnalyzer::write(char* outXMLFileName) 
{
    // Dispatch based on the file extension
    if (extension == ".vtu")
    {
      if (unsGridReader) {
        WriteAnXMLFile<vtkXMLUnstructuredGridWriter> (outXMLFileName, unsGridReader->GetOutput());
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
    else if (extension == ".vto")
    {
      //dataSet = ReadAnXMLFile<vtkXMLHyperOctreeReader> (xmlFileName);
    }
    else if (extension == ".vtk")
    {
      if (dataSetReader) {
        WriteAnXMLFile<vtkXMLUnstructuredGridWriter> (outXMLFileName, dataSetReader->GetOutput());
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

double* vtkAnalyzer::getPointCoords(int pntId)
{
   double* pntCoords;
   if (pntId < dataSet->GetNumberOfPoints()) {
     pntCoords = dataSet->GetPoint(pntId);
   } else {
     std::cerr << "Point ID is out of range!" << std::endl;
   }
   return pntCoords;
}

// returns coordinates of member points in cell with ID
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
vtkAnalyzer::getCellCenters(int& numComponent, int nDim)
{
  int num_cells = getNumberOfCells();
  std::vector<double> cellCenters(num_cells*nDim,0);
  for (int i = 0; i < num_cells; ++i) {
    std::vector<double *> cellCoords = getCellCoords(i, numComponent);
    for (int j = 0; j < numComponent ; ++j) {
      cellCenters[i*nDim] += cellCoords[j][0]/numComponent;
      cellCenters[i*nDim+1] += cellCoords[j][1]/numComponent;
      cellCenters[i*nDim+2] += cellCoords[j][2]/numComponent;
      delete cellCoords[j];
    }
  }
  return cellCenters;
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
                                  std::ostream& outputStream)
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
    outputStream << std::left << std::setw(10) << i << std::left << std::setw(16); 
    for (int j = 0; j < interpData.size(); ++j) {
      outputStream << interpData[j][i] << std::left << std::setw(16);
    }
    double G = interpData[0][i]*R*T*(1 - Mc/M);
    // if V not given, get from G and E
    if (poisson_dom_default == -1) {
      double V = (youngs_dom_default*G)/(3*(3*G - youngs_dom_default));
    // correcting for negative 0 output
      outputStream << std::left << std::setw(16) << youngs_dom_default 
                   << std::left << std::setw(16) << (V == 0.0 ? abs(V): V) << std::endl;
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

// no spheres and coord switch
void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default,
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::ostream& outputStream, bool writeCoord)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  }
  if (!writeCoord)
    writeInterpData(interpData, Mc, M, youngs_dom_default,
                    poisson_dom_default, T, R,
                    outputStream);
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
      outputStream << std::left << std::setw(10) << i 
                   << std::left << std::setw(16) << PlaneCellCenters[i*nDim] 
                   << std::left << std::setw(16) << PlaneCellCenters[i*nDim+1] 
                   << std::left << std::setw(16) << PlaneCellCenters[i*nDim+2] ;
      for (int j = 0; j < interpData.size(); ++j) {
        outputStream << std::left << std::setw(16) << interpData[j][i]
                     << std::left << std::setw(16);
      }  
      double G = interpData[0][i]*R*T*(1 - Mc/M);
      // if V not given, get from G and E
      if (poisson_dom_default == -1) {
        double V = (youngs_dom_default*G)/(3*(3*G - youngs_dom_default));
      // correcting for negative 0 output
        outputStream << std::left << std::setw(16) << youngs_dom_default 
                     << std::left << std::setw(16) << (V == 0.0 ? abs(V): V) << std::endl;
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

// has spheres, write_coords=0
void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default,
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::vector<sphere>& spheres,
                                  std::vector<string>& mat_sphere_names,
                                  std::vector<string>& material_names,
                                  std::vector<double>& youngs_inc_default,
                                  std::vector<double>& shear_inc_default,
                                  std::vector<double>& poisson_inc_default,
                                  std::ostream& outputStream)
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
    outputStream << std::left << std::setw(10) << i << std::left << std::setw(16); 
    bool in_sphere = false;
    int k;
    for (int j = 0; j < interpData.size(); ++j) {
      // checking if plane points are in sphere
      std::vector<double> point;
      point.push_back(PlaneCellCenters[i*nDim]);
      point.push_back(PlaneCellCenters[i*nDim+1]);
      point.push_back(PlaneCellCenters[i*nDim+2]);
      for (k = 0; k < spheres.size(); ++k) {
        if(in_sphere=spheres[k].in_sphere(point)) {
          break;
        }
      }
      if (in_sphere)
        outputStream << 0.0 << std::left << std::setw(16);
      else      
        outputStream << interpData[j][i] << std::left << std::setw(16);
    }
    if (in_sphere) {
      for (int i = 0; i < material_names.size(); ++i) {
          if (!material_names[i].compare(mat_sphere_names[k])) {
            // if V not given, get from G and E
            if (poisson_inc_default[i] == -1) {
              outputStream << std::left << std::setw(16) << youngs_inc_default[i]
                           << std::left << std::setw(16)
                           << (youngs_inc_default[i]*shear_inc_default[i])/
                              (3*(3*shear_inc_default[i] - youngs_inc_default[i]))
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
      double G = interpData[0][i]*R*T*(1 - Mc/M);
      // if V not given, get from G and E
      if (poisson_dom_default == -1) {
      double V = (youngs_dom_default*G)/(3*(3*G - youngs_dom_default));
      // correcting for negative 0 output
      outputStream << std::left << std::setw(16) << youngs_dom_default
                   << std::left << std::setw(16) << (V == 0.0 ? abs(V): V) << std::endl;
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

// has spheres and coord switch
void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default, 
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::vector<sphere>& spheres,
                                  std::vector<string>& mat_sphere_names,
                                  std::vector<string>& material_names,
                                  std::vector<double>& youngs_inc_default,
                                  std::vector<double>& shear_inc_default,
                                  std::vector<double>& poisson_inc_default,
                                  std::ostream& outputStream, bool writeCoord)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  } 
  //double R = .000008314;
  //double T = 300.0; 
  if (!writeCoord)
    writeInterpData(interpData, Mc, M, youngs_dom_default, 
                    poisson_dom_default, T, R,
                    PlaneCellCenters, nDim, spheres, 
                    mat_sphere_names, material_names, 
                    youngs_inc_default, shear_inc_default, 
                    poisson_inc_default, outputStream);
  else {

    outputStream << std::left << std::setw(10) << "id" 
                 << std::left << std::setw(16) << "X"
                 << std::left << std::setw(16) << "Y"
                 << std::left << std::setw(16) << "Z" 
                 << std::left << std::setw(16) << "rho"   
                 << std::left << std::setw(16) << "E"
                 << std::left << std::setw(16) << "V" << std::endl << std::endl;
    
    for (int i = 0; i < interpData[0].size(); ++i) {
      outputStream << std::left << std::setw(10) << i  
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
          if(in_sphere=spheres[k].in_sphere(point)) {
            break;
          }
        }
        // if point in plane is in sphere, crosslink = 0
        if (in_sphere) 
          outputStream << std::left << std::setw(16) << 0.0 
                       << std::left << std::setw(16);
        else      
          outputStream << std::left << std::setw(16) << interpData[j][i] 
                       << std::left << std::setw(16);
      } 
      if (in_sphere) {
        for (int i = 0; i < material_names.size(); ++i) { 
          if (!material_names[i].compare(mat_sphere_names[k])) {
            // if V not given, get from G and E
            if (poisson_inc_default[i] == -1) {
              outputStream << std::left << std::setw(16) << youngs_inc_default[i] 
                           << std::left << std::setw(16) 
                           << (youngs_inc_default[i]*shear_inc_default[i])/
                              (3*(3*shear_inc_default[i] - youngs_inc_default[i]))
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
        double G = interpData[0][i]*R*T*(1 - Mc/M);
        // if V not given, get from G and E
        if (poisson_dom_default == -1) {
          double V = (youngs_dom_default*G)/(3*(3*G - youngs_dom_default));
        // correcting for negative 0 output
          outputStream << std::left << std::setw(16) << youngs_dom_default 
                       << std::left << std::setw(16) << (V == 0.0 ? abs(V): V) << std::endl;
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



void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                                  double Mc, double M, double youngs_dom_default, 
                                  double poisson_dom_default, double T, double R,
                                  const std::vector<double>& PlaneCellCenters, int nDim,
                                  std::string filename, bool writeCoord)
{
  std::ofstream outputStream(filename.c_str());
  if(!outputStream.good()) {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
  writeInterpData(interpData, Mc, M, youngs_dom_default, 
                  poisson_dom_default, T, R, PlaneCellCenters, nDim, 
                  outputStream, writeCoord);
}                                 


void vtkAnalyzer::writeInterpData(const std::vector<std::vector<double>>& interpData,
                       double Mc, double M, double youngs_dom_default, 
                       double poisson_dom_default, double T, double R,
                       const std::vector<double>& PlaneCellCenters, int nDim,
                       std::vector<sphere>& spheres,
                       std::vector<string>& mat_sphere_names,
                       std::vector<string>& material_names,
                       std::vector<double>& youngs_inc_default,
                       std::vector<double>& shear_inc_default,
                       std::vector<double>& poisson_inc_default,
                       std::string filename, bool writeCoord)
{
  std::ofstream outputStream(filename.c_str());
  if(!outputStream.good()) {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
  writeInterpData(interpData, Mc, M, youngs_dom_default, poisson_dom_default, T, R,
                  PlaneCellCenters, nDim, spheres, mat_sphere_names, material_names, 
                  youngs_inc_default, shear_inc_default, poisson_inc_default, 
                  outputStream, writeCoord);
}                                 

void vtkAnalyzer::writeCSV(char* fname, std::vector<double> slnVec)
{
   ofstream csvFile;
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
  vtkDoubleArray* data;
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
        double* comps;
        // EDIT: causing indexing issues
        //       shouldn't resize if pushing back
        //  pntData[iTup].resize(numComponent);
        comps = da->GetTuple(iTup);
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
   vtkDoubleArray* data;
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


void vtkAnalyzer::setPointDataArray(const char* name, int numComponent, 
                                    std::vector<double> &pntArrData)
{
   if (!pointData)
    pointData = dataSet->GetPointData();
   if (pointData)
   { 
   vtkSmartPointer<vtkDoubleArray> da = 
    vtkSmartPointer<vtkDoubleArray>::New();
   da->SetName(name);
   da->SetNumberOfComponents(numComponent);
   for(int i=0; i<getNumberOfPoints(); i++)
      da->InsertNextTuple(&pntArrData[i]);
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
      da->InsertNextTuple(&cellArrData[i]);
   cellData->SetActiveScalars(name);
   cellData->SetScalars(da);
   }
}

//-----------------------------------------------------------------------------
