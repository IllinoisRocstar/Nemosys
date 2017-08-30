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
													 std::vector<std::vector<double>> volDataMat,
													 std::vector<double> PlaneCellCenters)
{
	std::vector<double> VolPointCoords = getAllPointCoords(nDim);
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
     	                         PlaneCellCenters, volData, interpData[j]);
	}
	delete VolPointInterp;
	return interpData;
}

// interpolate point data from 3D mesh in neighborhoods of
// cell centers of planar mesh to those centers for cases with spheres
// numTuple should = VolPointCoords.size()/ndim
std::vector<std::vector<double>>
vtkAnalyzer::getInterpData(int nDim, int num_neighbors, int numComponent, int numTuple,
                           std::vector<std::vector<double>> volDataMat,
                           std::vector<double> PlaneCellCenters,
													 std::vector<sphere> spheres)
{
  std::vector<double> VolPointCoords = getAllPointCoords(nDim);
  int num_vol_points = getNumberOfPoints(); 
  int num_interp_points = PlaneCellCenters.size()/nDim; 

  basicInterpolant* VolPointInterp = 
    new basicInterpolant(nDim, num_vol_points, num_neighbors, VolPointCoords);

  std::vector<std::vector<double>> interpData(numComponent);  
  for (int j = 0; j < numComponent; ++j) {   
    std::vector<double > volData(numTuple);
    for (int i = 0; i < numTuple; ++i) {
    	bool in_sphere = false;
			// checking if vol points are in sphere
			std::vector<double> point;
			point.push_back(VolPointCoords[i*nDim]);
			point.push_back(VolPointCoords[i*nDim+1]);
			point.push_back(VolPointCoords[i*nDim+2]);
			for (int k = 0; k < spheres.size(); ++k) {
				if (spheres[i].in_sphere(point)) {
					in_sphere = true;
					break;
				}
			}
			if (in_sphere)
				volData[i] = 0;
			else	
				volData[i] = volDataMat[i][j];
    	}   
    VolPointInterp->interpolate(num_interp_points, 
                               PlaneCellCenters, volData, interpData[j]);
  }
  delete VolPointInterp;
  return interpData;
}


void vtkAnalyzer::writeInterpData(std::vector<std::vector<double>> interpData,
																	double Mc, double M, double E,
																	std::ostream& outputStream)
{
	if (!outputStream.good()) {
		std::cout << "Output stream is bad" << std::endl;
		exit(1);
	}

	double R = .000008314;
	double T = 300.0;	
	for (int i = 0; i < interpData[0].size(); ++i) {
    outputStream << i << "\t"; 
    for (int j = 0; j < interpData.size(); ++j) {
      outputStream << interpData[j][i] << "\t";
    }
		double G = interpData[0][i]*R*T*(1 - Mc/M);
		double V = (E*G)/(3*(3*G - E));
		outputStream << E << "\t" << V << std::endl;	
  }   
}


void vtkAnalyzer::writeInterpData(std::vector<std::vector<double>> interpData,
																	double Mc, double M, double E, 
																	std::vector<double> PlaneCellCenters, int nDim,
																	std::vector<sphere> spheres,
																	std::ostream& outputStream)
{
	if (!outputStream.good()) {
		std::cout << "Output stream is bad" << std::endl;
		exit(1);
	}
  
	double R = .000008314;
	double T = 300.0;	
	for (int i = 0; i < interpData[0].size(); ++i) {
    outputStream << i << "\t"; 
    for (int j = 0; j < interpData.size(); ++j) {
			// checking if plane points are in sphere
			bool in_sphere = false;
			std::vector<double> point;
			point.push_back(PlaneCellCenters[i*nDim]);
			point.push_back(PlaneCellCenters[i*nDim+1]);
			point.push_back(PlaneCellCenters[i*nDim+2]);
			for (int k = 0; k < spheres.size(); ++k) {
				if(spheres[i].in_sphere(point)) {
					in_sphere=true;
					break;
				}
			}
			if (in_sphere)
				outputStream << 0 << "\t";
			else			
				outputStream << interpData[j][i] << "\t";
    }
		double G = interpData[0][i]*R*T*(1 - Mc/M);
		double V = (E*G)/(3*(3*G - E));
		outputStream << E << "\t" << V << std::endl;	
  }   
}
void vtkAnalyzer::writeInterpData(std::vector<std::vector<double>> interpData,
																	double Mc, double M, double E,
                                  std::vector<double> PlaneCellCenters, int nDim,
                                  std::vector<sphere> spheres,
                                  std::ostream& outputStream, bool writeCoord)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  }
	double R = .000008314;
	double T = 300.0;	
  if (!writeCoord)
		writeInterpData(interpData, Mc, M, E, PlaneCellCenters, nDim, spheres, outputStream);
	else {
  	for (int i = 0; i < interpData[0].size(); ++i) {
    	outputStream << i << "\t"
      	           << "(" << PlaneCellCenters[i*nDim] << ","
                          << PlaneCellCenters[i*nDim+1] << ","
         	                << PlaneCellCenters[i*nDim+2] << ")";
    	for (int j = 0; j < interpData.size(); ++j) {
      	// checking if plane points are in sphere
      	bool in_sphere = false;
      	std::vector<double> point;
      	point.push_back(PlaneCellCenters[i*nDim]);
      	point.push_back(PlaneCellCenters[i*nDim+1]);
      	point.push_back(PlaneCellCenters[i*nDim+2]);
      	for (int k = 0; k < spheres.size(); ++k) {
        	if(spheres[i].in_sphere(point)) {
          	in_sphere=true;
          	break;
        	}
      	}
      	if (in_sphere) 
        	outputStream << "\t" << 0 << "\t";
      	else      
        	outputStream << "\t" << interpData[j][i] << "\t";
    	}  
			double G = interpData[0][i]*R*T*(1 - Mc/M);
			double V = (E*G)/(3*(3*G - E));
			outputStream << E << "\t" << V << std::endl;	
  	}
	}
}

void vtkAnalyzer::writeInterpData(std::vector<std::vector<double>> interpData,
																	double Mc, double M, double E,
																	std::vector<double> PlaneCellCenters, int nDim,
																	std::ostream& outputStream)
{
  if (!outputStream.good()) {
    std::cout << "Output stream is bad" << std::endl;
    exit(1);
  }
  
	double R = .000008314;
	double T = 300.0;	
  for (int i = 0; i < interpData[0].size(); ++i) {
	  outputStream << i << "\t"
    						 << "(" << PlaneCellCenters[i*nDim] << ","
                          << PlaneCellCenters[i*nDim+1] << ","
                          << PlaneCellCenters[i*nDim+2] << ")";
    for (int j = 0; j < interpData.size(); ++j) {
      outputStream << "\t" << interpData[j][i] << "\t";
    }  
		double G = interpData[0][i]*R*T*(1 - Mc/M);
		double V = (E*G)/(3*(3*G - E));
		outputStream << E << "\t" << V << std::endl;	
  }   

}

void vtkAnalyzer::writeInterpData(std::vector<std::vector<double>> interpData, 
																	double Mc, double M, double E,
																  std::string filename)
{
	std::ofstream outputStream(filename.c_str());
	if (!outputStream.good()) {
		std::cout << "Output file stream is bad" << std::endl;
		exit(1);
	}
	writeInterpData(interpData, Mc, M, E, outputStream);
}

void vtkAnalyzer::writeInterpData(std::vector<std::vector<double>> interpData,
																	double Mc, double M, double E,
																	std::vector<double> PlaneCellCenters, int nDim,
																	std::string filename)
{
	std::ofstream outputStream(filename.c_str());
	if(!outputStream.good()) {
		std::cout << "Output file stream is bad" << std::endl;
		exit(1);
	}
	writeInterpData(interpData, Mc, M, E, PlaneCellCenters, nDim, outputStream);
}																	


void vtkAnalyzer::writeInterpData(std::vector<std::vector<double>> interpData,
																	double Mc, double M, double E,
																	std::vector<double> PlaneCellCenters, int nDim,
                                  std::vector<sphere> spheres,
                                  std::string filename, bool writeCoord)
{
	std::ofstream outputStream(filename.c_str());
	if(!outputStream.good()) {
		std::cout << "Output file stream is bad" << std::endl;
		exit(1);
	}
	writeInterpData(interpData, Mc, M, E, PlaneCellCenters, nDim, spheres, outputStream, writeCoord);
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
				//	pntData[iTup].resize(numComponent);
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
