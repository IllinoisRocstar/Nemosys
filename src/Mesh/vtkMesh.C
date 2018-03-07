#include<vtkMesh.H>

void vtkMesh::write()
{
  if (!dataSet)
  {
    std::cout << "No dataSet to write!" << std::endl;
    exit(1);
  }
   
  std::string extension = find_ext(filename);
 
  if (extension == ".vtp")
    writeVTFile<vtkXMLPolyDataWriter> (filename,dataSet);
  else if (extension == ".vts")
    writeVTFile<vtkXMLStructuredGridWriter> (filename, dataSet);
  else if (extension == ".vtr")
    writeVTFile<vtkXMLRectilinearGridWriter> (filename,dataSet);
  else if (extension == ".vti")
    writeVTFile<vtkXMLImageDataWriter> (filename, dataSet);
  else if (extension == ".vto")
    writeVTFile<vtkXMLHyperOctreeWriter> (filename,dataSet);
  else
    writeVTFile<vtkXMLUnstructuredGridWriter> (filename,dataSet);   // default is vtu 
} 

void vtkMesh::write(std::string fname)
{
  if (!dataSet)
  {
    std::cout << "No dataSet to write!" << std::endl;
    exit(1);
  }
   
  std::string extension = find_ext(fname);
 
  if (extension == ".vtp")
    writeVTFile<vtkXMLPolyDataWriter> (fname,dataSet);
  else if (extension == ".vts")
    writeVTFile<vtkXMLStructuredGridWriter> (fname, dataSet);
  else if (extension == ".vtr")
    writeVTFile<vtkXMLRectilinearGridWriter> (fname,dataSet);
  else if (extension == ".vti")
    writeVTFile<vtkXMLImageDataWriter> (fname, dataSet);
  else if (extension == ".vto")
    writeVTFile<vtkXMLHyperOctreeWriter> (fname,dataSet);
  else
    writeVTFile<vtkXMLUnstructuredGridWriter> (fname,dataSet);   // default is vtu 
  
}

vtkMesh::vtkMesh(const char* fname)
{
  std::string extension = vtksys::SystemTools::GetFilenameLastExtension(fname);
  // Dispatch based on the file extension
  if (extension == ".vtu")
  {
    dataSet.TakeReference(ReadAnXMLFile<vtkXMLUnstructuredGridReader> (fname));
  }
  else if (extension == ".vtp")
  {
    dataSet.TakeReference(ReadAnXMLFile<vtkXMLPolyDataReader> (fname));
  }
  else if (extension == ".vts")
  {
    dataSet.TakeReference(ReadAnXMLFile<vtkXMLStructuredGridReader> (fname));
  }
  else if (extension == ".vtr")
  {
    dataSet.TakeReference(ReadAnXMLFile<vtkXMLRectilinearGridReader> (fname));
  }
  else if (extension == ".vti")
  {
    dataSet.TakeReference(ReadAnXMLFile<vtkXMLImageDataReader> (fname));
  }
  else if (extension == ".vto")
  {
    dataSet.TakeReference(ReadAnXMLFile<vtkXMLHyperOctreeReader> (fname));
  }
  else if (extension == ".vtk")
  {
    //dataSet.TakeReference(ReadALegacyVTKFile(fname));
    dataSet = vtkDataSet::SafeDownCast(ReadALegacyVTKFile(fname));
  }
  else
  {
    std::cerr << "Unknown extension: " << extension << std::endl;
    exit(1);
  }

  if (!dataSet)
  {
    std::cout << "Error populating dataSet" << std::endl;
    exit(1);
  }

  std::string newname(fname);
  setFileName(newname);
  std::cout << "vtkMesh constructed" << std::endl;
  numCells = dataSet->GetNumberOfCells();
  numPoints = dataSet->GetNumberOfPoints();
}

vtkSmartPointer<vtkUnstructuredGrid> ReadALegacyVTKFile(const char* fileName)
{
  std::ifstream meshStream(fileName);
  if (!meshStream.good())
  {
    std::cerr << "Could not open file " << fileName << std::endl;
    exit(1);
  } 

  // declare points to be pushed into dataSet_tmp
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  // declare vector of cell ids to be pushed into dataSet_tmp
  std::vector<vtkSmartPointer<vtkIdList>> vtkCellIds;
  // declare dataSet_tmp which will be associated to output vtkMesh
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = vtkSmartPointer<vtkUnstructuredGrid>::New();
  int numPoints, numCells;
  double point[3];
  
  std::string line;
  while (getline(meshStream,line))
  {
    // assert that file contains an unstructured grid
    if (line.find("DATASET") != -1)
    {
      if (line.find("UNSTRUCTURED_GRID") == -1)
      {
        std::cerr << "Reading a " << line << " is not supported" << std::endl;
        exit(1);
      }
    }

    if (line.find("FIELD") != -1)
    {
      std::cout <<" in field" << std::endl;
      {
        std::stringstream ss(line);
        std::string dataname;
        int numarr;
        ss >> dataname >> numarr;
      }
      getline(meshStream,line);
      std::stringstream ss(line);
      std::string arrname, type;
      int numComponent, numTuple;
      ss >> arrname >> numComponent >> numTuple >> type;
      // TODO: extend this to support other types
      if (!type.compare("int"))
      {
        vtkSmartPointer<vtkIntArray> arr = vtkSmartPointer<vtkIntArray>::New();
        arr->SetName(&arrname[0u]);
        arr->SetNumberOfComponents(numComponent);
        arr->SetNumberOfTuples(numTuple);
        getline(meshStream,line);
        std::stringstream ss(line);
        int val;
        for (int i = 0; i < numTuple; ++i)
        {
          ss >> val;
          arr->InsertValue(i,val);
        }
        dataSet_tmp->GetFieldData()->AddArray(arr);
      }  

    }    
 
    if (line.find("POINTS") != - 1)
    {
      std::stringstream ss(line);
      std::string tmp;
      ss >> tmp >> numPoints;
      points->SetNumberOfPoints(numPoints);
      for (int i = 0; i < numPoints; ++i)
      {
        getline(meshStream,line);
        stringstream ss(line);
        ss >> point[0] >> point[1] >> point[2];
        points->SetPoint(i,point);
      }
      dataSet_tmp->SetPoints(points);
    }

    if (line.find("CELLS") != -1)
    {
      std::stringstream ss(line);
      std::string tmp;
      ss >> tmp >> numCells;
      dataSet_tmp->Allocate(numCells);
      vtkCellIds.resize(numCells);
      for (int i = 0; i < numCells; ++i)
      {
        getline(meshStream,line);
        std::stringstream ss(line);
        int numId;
        double id;
        ss >> numId;
        vtkCellIds[i] = vtkSmartPointer<vtkIdList>::New();
        vtkCellIds[i]->SetNumberOfIds(numId);
        for (int j = 0; j < numId; ++j)
        {
          ss >> id;
          vtkCellIds[i]->SetId(j,id);
        } 
      } 
    }
    if (line.find("CELL_TYPES") != -1)
    {
      for (int i = 0; i < numCells; ++i)
      {
        getline(meshStream,line);
        std::stringstream ss(line);
        int cellType;
        ss >> cellType;
        dataSet_tmp->InsertNextCell(cellType,vtkCellIds[i]);  
      }
    }
  
    if (line.find("CELL_DATA") != -1)
    {
      getline(meshStream,line);
      std::stringstream ss(line);
      std::string attribute, name, type;
      int numComponent = 1;
      ss >> attribute >> name >> type >> numComponent;
      // TODO: extend by iterating over permutations of attribute types
      if (!attribute.compare("SCALARS"))
      {
        numComponent = (numComponent > 4 || numComponent < 1 ? 1 : numComponent);
        vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
        arr->SetName(&name[0u]);
        arr->SetNumberOfComponents(numComponent);
        arr->SetNumberOfTuples(numCells);
        // skip lookup table line
        getline(meshStream,line);
        for (int i = 0; i < numCells; ++i)
        {
          getline(meshStream,line);
          std::stringstream ss(line);
          double vals[numComponent];
          for (int j = 0; j < numComponent; ++j)
          {
            ss >> vals[j];
          }
          arr->InsertTuple(i,vals); 
        }
        dataSet_tmp->GetCellData()->AddArray(arr);
        getline(meshStream,line);
        std::stringstream ss(line);
        ss >> attribute >> name >> type >> numComponent; 
        if (!attribute.compare("SCALARS"))
        {
          numComponent = (numComponent > 4 || numComponent < 1 ? 1 : numComponent);
          vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
          arr->SetName(&name[0u]);
          arr->SetNumberOfComponents(numComponent);
          arr->SetNumberOfTuples(numCells);
          // skip lookup table line
          getline(meshStream,line);
          for (int i = 0; i < numCells; ++i)
          {
            getline(meshStream,line);
            std::stringstream ss(line);
            double vals[numComponent];
            for (int j = 0; j < numComponent; ++j)
            {
              ss >> vals[j];
            }
            arr->InsertTuple(i,vals); 
          }
          dataSet_tmp->GetCellData()->AddArray(arr);
        }  
      } 
    }
    if (line.find("POINT_DATA")  != -1)
    {
      getline(meshStream,line); 
      std::stringstream ss(line);
      std::string attribute, name, type;
      ss >> attribute >> name >> type;
      if (!attribute.compare("VECTORS"))
      {
        vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
        arr->SetName(&name[0u]);
        arr->SetNumberOfComponents(3);
        arr->SetNumberOfTuples(numPoints);
        for (int i = 0; i < numPoints; ++i)
        {
          getline(meshStream, line);
          std::stringstream ss(line);
          double vec[3];
          for (int j = 0; j < 3; ++j)
          {
            ss >> vec[j];
          } 
          arr->InsertTuple(i,vec);
        }
        dataSet_tmp->GetPointData()->AddArray(arr);
        getline(meshStream,line); 
        std::stringstream ss(line);
        std::string attribute, name, type;
        ss >> attribute >> name >> type;
        if (!attribute.compare("VECTORS"))
        {
          vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
          arr->SetName(&name[0u]);
          arr->SetNumberOfComponents(3);
          arr->SetNumberOfTuples(numPoints);
          for (int i = 0; i < numPoints; ++i)
          {
            getline(meshStream, line);
            std::stringstream ss(line);
            double vec[3];
            for (int j = 0; j < 3; ++j)
            {
              ss >> vec[j];
            } 
            arr->InsertTuple(i,vec);
          }
          dataSet_tmp->GetPointData()->AddArray(arr);
          
        }
      } 
    }  

  }
  return dataSet_tmp;
}

// get point with id
std::vector<double> vtkMesh::getPoint(int id)
{
  double coords[3];
  dataSet->GetPoint(id, coords);
  std::vector<double> result(coords, coords+3);
  return result;
}

// get cell with id : returns point indices and respective coordinates
std::map<int, std::vector<double>> vtkMesh::getCell(int id)
{
  if (id < numCells) 
  {
    std::map<int,std::vector<double>> cell;
    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    point_ids = dataSet->GetCell(id)->GetPointIds();
    int num_ids = point_ids->GetNumberOfIds();
    for (int i = 0; i < num_ids; ++i) 
    {
      int pntId = point_ids->GetId(i);
      std::vector<double> coord = getPoint(pntId);
      cell.insert(std::pair<int,std::vector<double>> (pntId,coord)); 
    }

    return cell;
  }
  else {
    std::cerr << "Cell ID is out of range!" << std::endl;
    exit(1);
  }
}

std::vector<std::vector<double>> vtkMesh::getCellVec(int id)
{
  if (id < numCells) 
  {
    std::vector<std::vector<double>> cell;
    vtkSmartPointer<vtkIdList> point_ids = vtkSmartPointer<vtkIdList>::New();
    point_ids = dataSet->GetCell(id)->GetPointIds();
    int num_ids = point_ids->GetNumberOfIds();
    cell.resize(num_ids);
    for (int i = 0; i < num_ids; ++i) 
    {
      int pntId = point_ids->GetId(i);
      cell[i] = getPoint(pntId);
    }
    return cell;
  }
  else {
    std::cerr << "Cell ID is out of range!" << std::endl;
    exit(1);
  }

}

void vtkMesh::report()
{
  if (!dataSet)
  {
    std::cout << "dataSet has not been populated!" << std::endl;
    exit(1);
  }

  typedef std::map<int,int> CellContainer;
  // Generate a report
  std::cout << "Processing the dataset generated from " << filename << std::endl
     << " dataSet contains a " 
     << dataSet->GetClassName()
     <<  " that has " << numCells << " cells"
     << " and " << numPoints << " points." << std::endl;

  CellContainer cellMap;
  for (int i = 0; i < numCells; i++)
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
      std::cout << "\tArray " << i << " is named "
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
    std::cout << " contains cell data with " << cd->GetNumberOfArrays()
              << " arrays." << std::endl;
    for (int i = 0; i < cd->GetNumberOfArrays(); i++)
    {
      std::cout << "\tArray " << i << " is named "
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

// get diameter of circumsphere of each cell
std::vector<double> vtkMesh::getCellLengths() 
{
  std::vector<double> result;
  result.resize(getNumberOfCells());
  for (int i = 0; i < getNumberOfCells(); ++i)
  {
    result[i] = std::sqrt(dataSet->GetCell(i)->GetLength2());
  } 
  return result;
}

// get center of a cell
std::vector<double> vtkMesh::getCellCenter(int cellID)
{
  std::vector<double> center(3);
  std::vector<std::vector<double>> cell = getCellVec(cellID); 
 
  for (int i = 0; i < cell.size(); ++i)
  {
    center = center + cell[i];
  }
  return (1./cell.size())*center;
}


// set point data (numComponets per point determined by dim of data[0] 
void vtkMesh::setPointDataArray(const char* name, const std::vector<std::vector<double>>& data)
{
  vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
  da->SetName(name);
  da->SetNumberOfComponents(data[0].size());
  for(int i=0; i < numPoints; i++)
    da->InsertNextTuple(data[i].data());
  dataSet->GetPointData()->AddArray(da);
  //dataSet->GetPointData()->SetActiveScalars(name);
  //dataSet->GetPointData()->SetScalars(da);
}

// set cell data (numComponents per cell determined by dim of data[0])
void vtkMesh::setCellDataArray(const char* name, const std::vector<std::vector<double>>& data)
{
  vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
  da->SetName(name);
  da->SetNumberOfComponents(data[0].size());
  for(int i=0; i < numCells; i++)
    da->InsertNextTuple(data[i].data());
  dataSet->GetCellData()->AddArray(da);
  //dataSet->GetCellData()->SetActiveScalars(name);
  //dataSet->GetCellData()->SetScalars(da);
}

void vtkMesh::setCellDataArray(const char* name, const std::vector<double>& data)
{   
  vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
  da->SetName(name);
  da->SetNumberOfComponents(1);
  for(int i=0; i < numCells; i++)
    da->InsertNextTuple1(data[i]);
  dataSet->GetCellData()->AddArray(da);
  //dataSet->GetCellData()->SetActiveScalars(name);
  //dataSet->GetCellData()->SetScalars(da);
}

// remove point data with given id from dataSet if it exists
void vtkMesh::unsetPointDataArray(int arrayID)
{
  dataSet->GetPointData()->RemoveArray(arrayID); 
}

void vtkMesh::unsetPointDataArray(const char* name)
{
  dataSet->GetPointData()->RemoveArray(name); 
}
// remove cell data with given id from dataSet if it exists
void vtkMesh::unsetCellDataArray(int arrayID)
{
  dataSet->GetCellData()->RemoveArray(arrayID);
}

void vtkMesh::unsetCellDataArray(const char* name)
{
  dataSet->GetCellData()->RemoveArray(name);
}

// remove field data with given id from dataSet
void vtkMesh::unsetFieldDataArray(const char* name)
{
  dataSet->GetFieldData()->RemoveArray(name);
}

