#include "vtkAnalyzer.H"

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
      //dataSet = ReadAnXMLFile<vtkXMLImageDataReader> (xmlFileName);
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
   if (pntId <= dataSet->GetNumberOfPoints()) {
     pntCoords = dataSet->GetPoint(pntId);
   } else {
     std::cerr << "Point ID is out of range!" << std::endl;
   }
   return pntCoords;
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
   if (pointData)
   {
   vtkDataArray* da = pointData->GetArray(id);
   if (da) {
      numComponent = da->GetNumberOfComponents();
      numTuple = da->GetNumberOfTuples();
      pntData.resize(numTuple);
      for (int iTup=0; iTup<numTuple; iTup++){
	  double* comps;
	  pntData[iTup].resize(numComponent);
	  comps = da->GetTuple(iTup);
	  for (int iComp=0; iComp<numComponent; iComp++){
	    pntData[iTup].push_back(comps[iComp]);
	  }
      }
      return(numComponent*numTuple);
   } else {
      return(0);
   }
   } else {
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
