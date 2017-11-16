#include<meshBase.H>
#include<vtkMesh.H>

template<class TReader> vtkDataSet* ReadAnXMLFile(const char* fileName)
{
  vtkSmartPointer<TReader> reader =
    vtkSmartPointer<TReader>::New();
  reader->SetFileName(fileName);
  reader->Update();
  reader->GetOutput()->Register(reader);
  return vtkDataSet::SafeDownCast(reader->GetOutput());
}

vtkMesh::vtkMesh(const char* fname)
{
  std::string extension = vtksys::SystemTools::GetFilenameLastExtension(fname);
  // Dispatch based on the file extension
  if (extension == ".vtu")
  {
    dataSet = ReadAnXMLFile<vtkXMLUnstructuredGridReader> (fname);
  }
  else if (extension == ".vtp")
  {
    dataSet = ReadAnXMLFile<vtkXMLPolyDataReader> (fname);
  }
  else if (extension == ".vts")
  {
    dataSet = ReadAnXMLFile<vtkXMLStructuredGridReader> (fname);
  }
  else if (extension == ".vtr")
  {
    dataSet = ReadAnXMLFile<vtkXMLRectilinearGridReader> (fname);
  }
  else if (extension == ".vti")
  {
    dataSet = ReadAnXMLFile<vtkXMLImageDataReader> (fname);
  }
  else if (extension == ".vto")
  {
    dataSet = ReadAnXMLFile<vtkXMLHyperOctreeReader> (fname);
  }
  else if (extension == ".vtk")
  {
    dataSet = ReadAnXMLFile<vtkDataSetReader> (fname);
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

  std::cout << "vtkMesh constructed" << std::endl;
  numCells = dataSet->GetNumberOfCells();
  numPoints = dataSet->GetNumberOfPoints();
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


