#include <vtkMesh.H>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkTriangleFilter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtksys/SystemTools.hxx>
#include <vtkCellTypes.h>
#include <vtkXMLWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkSTLWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkSTLReader.h>
#include <vtkExtractEdges.h>
#include <vtkGenericCell.h>
#include <vtkStructuredGrid.h>
#include <vtkCellIterator.h>
#include <vtkRectilinearGrid.h>
#include <vtkImageData.h>
#include <AuxiliaryFunctions.H>

using nemAux::operator*; // for vector multiplication.
using nemAux::operator+; // for vector addition.

void vtkMesh::write(const std::string &fname) const
{
  if (!dataSet)
  {
    std::cout << "No dataSet to write!" << std::endl;
    exit(1);
  }
   
  std::string extension = nemAux::find_ext(fname);

  if (extension == ".vtp")
    writeVTFile<vtkXMLPolyDataWriter> (fname,dataSet);
  else if (extension == ".vts")
    writeVTFile<vtkXMLStructuredGridWriter> (fname, dataSet);
  else if (extension == ".vtr")
    writeVTFile<vtkXMLRectilinearGridWriter> (fname,dataSet);
  else if (extension == ".vti")
    writeVTFile<vtkXMLImageDataWriter> (fname, dataSet);
  else if (extension == ".stl")
  {
    writeVTFile<vtkSTLWriter> (fname, dataSet); // ascii stl
  }
  else if (extension == ".vtk")
  {
    writeVTFile<vtkUnstructuredGridWriter> (fname, dataSet); // legacy vtk writer
  }
  else {
    std::string fname_tmp = nemAux::trim_fname(fname, ".vtu");
    writeVTFile<vtkXMLUnstructuredGridWriter>(fname_tmp, dataSet);   // default is vtu
  }
}

vtkMesh::vtkMesh(vtkSmartPointer<vtkDataSet> dataSet_tmp,
                 const std::string &fname)
{
  if (dataSet_tmp) {
    dataSet = dataSet_tmp;
    filename = fname;
    numCells = dataSet->GetNumberOfCells();
    numPoints = dataSet->GetNumberOfPoints();
    std::cout << "vtkMesh constructed" << std::endl;
  } else {
    std::cerr << "Nothing to copy!" << std::endl;
    exit(1);
  }
}

vtkMesh::vtkMesh(const std::vector<double>& xCrds,
                 const std::vector<double>& yCrds,
                 const std::vector<double>& zCrds,
                 const std::vector<int>& elemConn,
                 const int cellType,
                 const std::string &newname)
{
  if (!(xCrds.size() == yCrds.size() && xCrds.size() == zCrds.size()))
  {
    std::cerr << "Length of coordinate arrays must match!" << std::endl;
    exit(1);
  }
  // point to be pushed into dataSet
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New(); 
  // declare vtk dataset
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp 
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  numPoints = xCrds.size();
  // allocate size for vtk point container 
  points->SetNumberOfPoints(numPoints);
  for (int i = 0; i < numPoints; ++i)
  {
    points->SetPoint(i, xCrds[i], yCrds[i], zCrds[i]); 
  }
  // add points to vtk mesh data structure 
  dataSet_tmp->SetPoints(points);  
  switch(cellType)
  {
    case VTK_TETRA:
    {
      numCells = elemConn.size()/4;
      dataSet_tmp->Allocate(numCells);
      for (int i = 0; i < numCells; ++i)
      {
        vtkSmartPointer<vtkIdList> elmConn = vtkSmartPointer<vtkIdList>::New();
        elmConn->SetNumberOfIds(4);
        for (int j = 0; j < 4; ++j)
        {
          elmConn->SetId(j,elemConn[i*4 + j]);
        }
        dataSet_tmp->InsertNextCell(VTK_TETRA, elmConn);
      }
      break;
    }
    case VTK_TRIANGLE:
    {
      numCells = elemConn.size()/3;
      dataSet_tmp->Allocate(numCells);
      for (int i = 0; i < numCells; ++i)
      {
        vtkSmartPointer<vtkIdList> elmConn = vtkSmartPointer<vtkIdList>::New();
        elmConn->SetNumberOfIds(3);
        for (int j = 0; j < 3; ++j)
        {
          elmConn->SetId(j,elemConn[i*3 + j]);
        }
        dataSet_tmp->InsertNextCell(VTK_TRIANGLE, elmConn);
      }
      break;
    }
    default:
    {
      std::cout << "Unknown element type " << cellType << std::endl;
      exit(1);
    }
  }
  filename = newname;  
  dataSet = dataSet_tmp; 
  std::cout << "vtkMesh constructed" << std::endl;
}

vtkMesh::vtkMesh(const char* fname)
{
  bool degenerate(0);
  std::string extension = vtksys::SystemTools::GetFilenameLastExtension(fname);
  // Dispatch based on the file extension
  if (extension == ".vtu")
  {
    dataSet.TakeReference(ReadAnXMLOrSTLFile<vtkXMLUnstructuredGridReader> (fname));
  }
  else if (extension == ".vtp")
  {
    dataSet.TakeReference(ReadAnXMLOrSTLFile<vtkXMLPolyDataReader> (fname));
  }
  else if (extension == ".vts")
  {
    dataSet.TakeReference(ReadAnXMLOrSTLFile<vtkXMLStructuredGridReader> (fname));
  }
  else if (extension == ".vtr")
  {
    dataSet.TakeReference(ReadAnXMLOrSTLFile<vtkXMLRectilinearGridReader> (fname));
  }
  else if (extension == ".vti")
  {
    dataSet.TakeReference(ReadAnXMLOrSTLFile<vtkXMLImageDataReader> (fname));
  }
  else if (extension == ".stl")
  {
    dataSet.TakeReference(ReadAnXMLOrSTLFile<vtkSTLReader> (fname));
  }
  else if (extension == ".vtk")
  {
    // if vtk is produced by MFEM, it's probably degenerate (i.e. point duplicity)
    // in this case, we use a different reader to correct duplicity and transfer 
    // the data attributes read by the regular legacy reader using the transfer
    // methods in meshBase
    std::ifstream meshStream(fname);
    if (!meshStream.good())
    {
      std::cout << "Error opening file " << fname << std::endl;
      exit(1);
    } 
    std::string line;
    getline(meshStream,line);
    getline(meshStream,line);
    meshStream.close();
    if (line.find("MFEM") != -1)
    {
      degenerate = 1;
      dataSet = vtkDataSet::SafeDownCast(ReadDegenerateVTKFile(fname));
      std::string newname(fname);
      setFileName(newname);
      numCells = dataSet->GetNumberOfCells();
      numPoints = dataSet->GetNumberOfPoints();
      vtkMesh vtkMesh_tmp(vtkDataSet::SafeDownCast(ReadALegacyVTKFile(fname)),fname);
      vtkMesh_tmp.transfer(this,"Consistent Interpolation");     
      std::cout << "vtkMesh constructed" << std::endl;
    }
    else
    {
      dataSet = vtkDataSet::SafeDownCast(ReadALegacyVTKFile(fname));
    }  

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

  if (!degenerate)
  {
    std::string newname(fname);
    setFileName(newname);
    std::cout << "vtkMesh constructed" << std::endl;
    numCells = dataSet->GetNumberOfCells();
    numPoints = dataSet->GetNumberOfPoints();
  }
}

vtkMesh::vtkMesh(const char* fname1, const char* fname2)
{
  std::string ext_in = vtksys::SystemTools::GetFilenameLastExtension(fname1);
  std::string ext_out = vtksys::SystemTools::GetFilenameLastExtension(fname2);

  // sanity check 
  if ( !(ext_in == ".vtu" && ext_out == ".stl") )
  {
    std::cerr << "vtkMesh: Currently only support conversion between vtu and stl.\n";
    exit(1);
  }

  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(fname1);
  reader->Update();
  reader->GetOutput()->Register(reader);

  // obtaining dataset
  dataSet.TakeReference(vtkDataSet::SafeDownCast(reader->GetOutput()));
  if (!dataSet)
  {
    std::cout << "Error populating dataSet" << std::endl;
    exit(1);
  }
  std::string newname(fname1);
  setFileName(newname);
  std::cout << "vtkMesh constructed" << std::endl;
  numCells = dataSet->GetNumberOfCells();
  numPoints = dataSet->GetNumberOfPoints();

  // skinn to the surface
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfFilt =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfFilt->SetInputConnection(reader->GetOutputPort());
  surfFilt->Update();

  // triangulate the surface
  vtkSmartPointer<vtkTriangleFilter> triFilt =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triFilt->SetInputConnection(surfFilt->GetOutputPort());
  triFilt->Update();

  // write to stl file
  vtkSmartPointer<vtkSTLWriter> stlWriter =
    vtkSmartPointer<vtkSTLWriter>::New();
  stlWriter->SetFileName(fname2);
  stlWriter->SetInputConnection(triFilt->GetOutputPort());
  stlWriter->Write();
}

bool readLegacyVTKHeader(const std::string& line)
{
  if (line.find("DATASET") != -1)
  {
    if (line.find("UNSTRUCTURED_GRID") == std::string::npos)
    {
      std::cerr << "Reading a " << line << " is not supported" << std::endl;
      exit(1);
    }
    return 1;
  }
  return 0;
}

bool readLegacyVTKFieldData(std::istream& meshStream, std::string& line, 
                            vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp)
{
  if (line.find("FIELD") != -1)
  {
    if (getline(meshStream,line))
    {
      std::string arrname, type;
      int numComponent, numTuple;
      std::stringstream ss;
      ss.str(line);
      ss >> arrname >> numComponent >> numTuple >> type;
      vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
      arr->SetName(&arrname[0u]);
      arr->SetNumberOfComponents(numComponent);
      arr->SetNumberOfTuples(numTuple);
      getline(meshStream,line);
      ss.str(line);
      double val;
      for (int i = 0; i < numTuple; ++i)
      {
        ss >> val;
        arr->InsertTuple(i,&val);
      }
      dataSet_tmp->GetFieldData()->AddArray(arr);
      return 1;
    }
  }
  return 0;
}

bool readLegacyVTKPoints(std::istream& meshStream, std::string& line, int& numPoints, 
                         vtkSmartPointer<vtkPoints> points,
                         vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp)
{
  if (line.find("POINTS") != -1)
  {
    std::istringstream ss(line);
    std::string tmp;
    ss >> tmp >> numPoints;
    points->SetNumberOfPoints(numPoints);
    double point[3];
    int i = 0;
    while (getline(meshStream,line) && i < numPoints)
    {
      if (!line.empty())
      {
        std::istringstream ss(line);
        ss >> point[0] >> point[1] >> point[2];
        points->SetPoint(i,point);
        ++i;
      }
    }
    dataSet_tmp->SetPoints(points);
    return 1;
  }
  return 0;
}

bool readLegacyVTKCells(std::istream& meshStream, std::string& line, int& numCells, 
                        std::vector<vtkSmartPointer<vtkIdList>>& vtkCellIds,
                        vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp)
{
  std::cout << line << std::endl;
  if (line.find("CELLS") != -1)
  {
    std::istringstream ss(line);
    std::string tmp;
    ss >> tmp >> numCells;
    vtkCellIds.resize(numCells);
    // get cells
    int i = 0;
    while (i < numCells && getline(meshStream,line))
    {
      if (!line.empty())
      {
        std::istringstream ss(line);
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
        ++i;
      } 
    }
    // get cell types
    while (getline(meshStream,line)) 
    {
      if (line.find("CELL_TYPES") != -1)
      {
        dataSet_tmp->Allocate(numCells);
        i = 0; 
        while (i < numCells && getline(meshStream,line))
        {
          if (!line.empty())
          { 
            std::istringstream ss(line);
            int cellType;
            ss >> cellType;
            dataSet_tmp->InsertNextCell(cellType,vtkCellIds[i]);  
            ++i;
          }
        }
        break; 
      }
    }
    return 1;
  }
  return 0;
}

bool readLegacyVTKData(std::ifstream& meshStream, std::string& line, const int numTuple, 
                       const bool pointOrCell, bool& hasPointOrCell,
                       vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp)
{
  // determines whether more than one point or cell array in dataSet
  bool inProgress(1); 
  // while there is a cell or point array to be read
  while (inProgress)
  {
    int numComponent(0);
    std::string line, attribute, name, type;
    vtkSmartPointer<vtkDoubleArray> arr = vtkSmartPointer<vtkDoubleArray>::New();
    while(getline(meshStream,line))
    {
      if (!line.empty())
      {
        if (line.find("LOOKUP_TABLE default") != std::string::npos)
        {
          break;
        }
        
        size_t hasScalar = line.find("SCALARS");
        size_t hasVector = line.find("VECTORS");
        size_t hasNormal = line.find("NORMALS");
        if (hasScalar == std::string::npos &&
            hasVector == std::string::npos &&
            hasNormal == std::string::npos)
        {
          std::cerr << "attribute type in " << line << " not supported" << std::endl;
          exit(1);
        }
        else
        { 
          std::istringstream ss(line);
          ss >> attribute >> name >> type >> numComponent;
          if (hasScalar != std::string::npos)
          { 
            numComponent = (numComponent > 4 || numComponent < 1) ? 1 : numComponent;
          }
          else if (hasVector != std::string::npos || hasNormal != std::string::npos)
          { 
            numComponent = 3;
          }
          arr->SetName(&name[0u]);
          arr->SetNumberOfComponents(numComponent);
          arr->SetNumberOfTuples(numTuple);
          //std::cout << arr->GetName() << " " << arr->GetNumberOfComponents() << " " 
          //          << arr->GetNumberOfTuples() << std::endl;
          break;
        }  
      }
    }

    int i = 0;
    while (i < numTuple && getline(meshStream,line))
    {
      if (!line.empty() && line.find("LOOKUP") == std::string::npos)
      {
        std::istringstream ss(line); 
        double vals[numComponent];
        for (int j = 0; j < numComponent; ++j)
        {
          ss >> vals[j];
        }
        arr->InsertTuple(i,vals); 
        ++i;
      }
    }
 
    if (pointOrCell)
      dataSet_tmp->GetPointData()->AddArray(arr);
    else
      dataSet_tmp->GetCellData()->AddArray(arr);

    bool moreData(0);
    int len(0);
    while (getline(meshStream,line))
    {
      len += line.length();
      // if point data encountered while looking for more cell data, we know we read all cell data
      if (line.find("POINT_DATA") != std::string::npos && !pointOrCell)
      {
        len -= line.length();
        hasPointOrCell = 1;
        break;
      }
      // if cell data encountered while looking for more point data, we know we read all point data
      if (line.find("CELL_DATA") != std::string::npos && pointOrCell)
      {
        hasPointOrCell = 1;
        len -= line.length();
        break;
      }
      // if the previous two conditions are not met and this one is, 
      // we know we have more than one point/cell data array
      if (line.find("SCALARS") != std::string::npos ||
          line.find("VECTORS") != std::string::npos ||
          line.find("NORMALS") != std::string::npos)
      {
        moreData = 1;
        break;   
      }
    }

    int curr = meshStream.tellg();
    meshStream.seekg(curr - len - 1); 

    if (!moreData)
      inProgress = 0;
  }
  return 1;
}

vtkSmartPointer<vtkUnstructuredGrid> ReadALegacyVTKFile(const char* fileName)
{
  std::ifstream meshStream;
  meshStream.open(fileName, std::ifstream::in);
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
  bool readHeader(0), readPoints(0), readCells(0), readFieldData(0), readCellData(0), readPointData(0);
  bool cellDataFirst(0), pointDataFirst(0), hasPointData(0), hasCellData(0); 
  std::string line;
  
  while(getline(meshStream, line))
  {
    if (!readHeader)
    {
      if (readHeader = readLegacyVTKHeader(line))
      { 
        std::cout << "read header" << std::endl;
      }
    }

    if (!readFieldData) 
    {
      if (readFieldData = readLegacyVTKFieldData(meshStream, line, dataSet_tmp))
      {
        std::cout << "read field data" << std::endl;
      }
    }    
 
    if (!readPoints) 
    {
      if (readPoints = readLegacyVTKPoints(meshStream, line, numPoints, points, dataSet_tmp))
      {
        std::cout << "read points" << std::endl;  
      }
    }

    if (!readCells) 
    {
      if (readCells = readLegacyVTKCells(meshStream, line, numCells, vtkCellIds, dataSet_tmp))
      {
        std::cout << "read cells" << std::endl;
      }
    }

    // call functions based on which data appears first in file
    if (!readCellData && !readPointData)
    {
      cellDataFirst 
        = (line.find("CELL_DATA") != std::string::npos && !pointDataFirst) ? 1 : 0;
      pointDataFirst  
        = (line.find("POINT_DATA") != std::string::npos && !cellDataFirst) ? 1: 0;
    }
    
    // read cell data then point data if it exists
    if (cellDataFirst)
    { 
      if (!readCellData) 
      {                      // boolean telling us if there is point data ------|
        if (readCellData = readLegacyVTKData(meshStream, line, numCells, 0, hasPointData, dataSet_tmp))
        {
          std::cout << "read cell data" << std::endl;
        }
      }
      
      if (!readPointData && hasPointData)
      {
        if(readPointData = readLegacyVTKData(meshStream, line, numPoints, 1, hasCellData, dataSet_tmp))
        {
          std::cout << "read point data" << std::endl;
        }
      }
    }
    // read point data and then cell data if it exists
    else if (pointDataFirst)
    {
      if (!readPointData)
      {
        if(readPointData = readLegacyVTKData(meshStream, line, numPoints, 1, hasCellData, dataSet_tmp))
        {
          std::cout << "read point data" << std::endl;
        }
      }
      
      if (!readCellData && hasCellData) 
      {
        if (readCellData = readLegacyVTKData(meshStream, line, numCells, 0, hasPointData, dataSet_tmp))
        {
          std::cout << "read cell data" << std::endl;
        }
      }
    }
  }
  return dataSet_tmp;
}

vtkSmartPointer<vtkUnstructuredGrid> ReadDegenerateVTKFile(const char* fileName)
{
  std::ifstream meshStream(fileName);
  if (!meshStream.good())
  {
    std::cerr << "Error opening file " << fileName << std::endl;
    exit(1);
  }
  std::string line;
  std::map<std::vector<double>, int> point_map;
  std::map<int, int> duplicate_point_map;
  std::pair<std::map<std::vector<double>,int>::iterator, bool> point_ret;
  int numPoints;
  int numCells;
  int cellListSize;
  std::vector<vtkSmartPointer<vtkIdList>> vtkCellIds;
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = vtkSmartPointer<vtkUnstructuredGrid>::New();
  while (getline(meshStream, line))
  {
    if (line.find("POINTS") != -1)
    {
      std::istringstream ss(line);
      std::string tmp;
      ss >> tmp >> numPoints;
      std::vector<double> point(3);
      int pntId = 0;
      int realPntId = 0;
      while (getline(meshStream,line) && pntId < numPoints)
      {
        if (!line.empty())
        {
          std::istringstream ss(line);
          ss >> point[0] >> point[1] >> point[2];
          point_ret = point_map.insert(std::pair<std::vector<double>,int> (point,realPntId));    
          if (point_ret.second)
          {
            ++realPntId;        
          }
          duplicate_point_map.insert(std::pair<int,int> (pntId,point_ret.first->second));
          ++pntId;
        }
      }
    }

    if (line.find("CELLS") != -1)
    {
      std::istringstream ss(line);
      std::string tmp;
      ss >> tmp >> numCells >> cellListSize;
      int cellId = 0;
      int realCellId = 0;
      vtkCellIds.resize(numCells);
      while (cellId < numCells && getline(meshStream,line))
      {
        if (!line.empty())
        {
          std::istringstream ss(line);
          int numId;
          double id;
          ss >> numId;
          vtkCellIds[cellId] = vtkSmartPointer<vtkIdList>::New();
          vtkCellIds[cellId]->SetNumberOfIds(numId);
          //cells[cellId].resize(numId);
          for (int j = 0; j < numId; ++j)
          {
            ss >> id;
            vtkCellIds[cellId]->SetId(j, duplicate_point_map[id]);
            //cells[cellId][j] = duplicate_point_map[id];
          }
          ++cellId;
        } 
      }
      // get cell types
      while (getline(meshStream,line)) 
      {
        if (line.find("CELL_TYPES") != -1)
        {
          dataSet_tmp->Allocate(numCells);
          cellId = 0;
          //cellTypes.resize(numCells); 
          while (cellId < numCells && getline(meshStream,line))
          {
            if (!line.empty())
            { 
              std::istringstream ss(line);
              int cellType;
              ss >> cellType;
              //cellTypes[cellId] = cellType;
              dataSet_tmp->InsertNextCell(cellType, vtkCellIds[cellId]);
              ++cellId;
            }
          }
          break; 
        }
      }
    }
  }
   
  std::multimap<int, std::vector<double>> flipped = nemAux::flip_map(point_map);
  std::multimap<int,std::vector<double>>::iterator flip_it = flipped.begin();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(point_map.size());
  int pnt = 0;
  while (flip_it != flipped.end())
  {
    points->SetPoint(pnt, flip_it->second.data());
    ++flip_it;
    ++pnt;
  }
  dataSet_tmp->SetPoints(points); 
  return dataSet_tmp;
}

// get point with id
std::vector<double> vtkMesh::getPoint(int id) const
{
  double coords[3];
  dataSet->GetPoint(id, coords);
  std::vector<double> result(coords, coords+3);
  return result;
}

// get 3 vectors with x,y and z coords
std::vector<std::vector<double>> vtkMesh::getVertCrds() const
{
  std::vector<std::vector<double>> comp_crds(3);
  for (int i = 0; i < 3; ++i)
  {
    comp_crds[i].resize(numPoints);
  }
  double coords[3];
  for (int i = 0; i < numPoints; ++i)
  {
    dataSet->GetPoint(i,coords);
    comp_crds[0][i] = coords[0];
    comp_crds[1][i] = coords[1];
    comp_crds[2][i] = coords[2];
  }
  return comp_crds;
}

// get cell with id : returns point indices and respective coordinates
std::map<int, std::vector<double>> vtkMesh::getCell(int id) const
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

std::vector<std::vector<double>> vtkMesh::getCellVec(int id) const
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


void vtkMesh::inspectEdges(const std::string& ofname) const
{
  std::ofstream outputStream(ofname);
  if (!outputStream.good())
  {
    std::cerr << "error opening " << ofname << std::endl;
    exit(1);
  }

  vtkSmartPointer<vtkExtractEdges> extractEdges 
    = vtkSmartPointer<vtkExtractEdges>::New();
  extractEdges->SetInputData(dataSet);
  extractEdges->Update();
  
  vtkSmartPointer<vtkGenericCell> genCell = vtkSmartPointer<vtkGenericCell>::New();
  for (int i = 0; i < extractEdges->GetOutput()->GetNumberOfCells(); ++i)
  {
    extractEdges->GetOutput()->GetCell(i, genCell);
    vtkPoints* points = genCell->GetPoints();
    double p1[3], p2[3];
    points->GetPoint(0,p1);
    points->GetPoint(1,p2);
    double len = sqrt(pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2));
    outputStream << len << std::endl;
  } 
}

std::vector<int> vtkMesh::getConnectivities() const
{
  std::vector<int> connectivities;
  vtkSmartPointer<vtkCellIterator> it 
    = vtkSmartPointer<vtkCellIterator>::Take(dataSet->NewCellIterator()); 
  for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell())
  {
    vtkSmartPointer<vtkIdList> pointIds = it->GetPointIds();
    for (int i = 0; i < pointIds->GetNumberOfIds(); ++i)
    {
      connectivities.push_back(pointIds->GetId(i));
    } 
  }
  return connectivities;
}

void vtkMesh::report() const
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

vtkSmartPointer<vtkDataSet> vtkMesh::extractSurface()
{
  // extract surface polygons
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfFilt =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfFilt->SetInputData(dataSet);
  surfFilt->Update();

  // triangulate the surface
  vtkSmartPointer<vtkTriangleFilter> triFilt =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triFilt->SetInputData(surfFilt->GetOutput());
  triFilt->Update();

  return triFilt->GetOutput();
}


// get diameter of circumsphere of each cell
std::vector<double> vtkMesh::getCellLengths() const
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
std::vector<double> vtkMesh::getCellCenter(int cellID) const
{
  std::vector<double> center(3,0.0);
  std::vector<std::vector<double>> cell = getCellVec(cellID);

  for (int i = 0; i < cell.size(); ++i)
  {
    center = center + cell[i];
  }
  return (1./cell.size())*center;
}

// returns the cell type
int vtkMesh::getCellType() const
{
  return dataSet->GetCellType(0);
}



// set point data  
void vtkMesh::setPointDataArray(const char* name, const std::vector<double>& data)
{
  vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
  da->SetName(name);
  da->SetNumberOfComponents(1);
  for(int i=0; i < numPoints; i++)
    da->InsertNextTuple1(data[i]);
  dataSet->GetPointData()->AddArray(da);
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

void vtkMesh::getPointDataArray(const std::string& name, std::vector<double>& data)
{
  int idx;  
  vtkSmartPointer<vtkDataArray> pd 
    = dataSet->GetPointData()->GetArray(&name[0u],idx);
  if (idx != -1)
  {
    if (pd->GetNumberOfComponents() > 1)
    {
      std::cerr << __func__ << " is only suitable for scalar data, i.e. 1 component\n";
      exit(1);
    }
    data.resize(pd->GetNumberOfTuples());
    double x[1];
    for (int i = 0; i < pd->GetNumberOfTuples(); ++i)
    {
      pd->GetTuple(i,x);
      data[i] = x[0];
    }
  }
  else
  {
    std::cerr << "could not find data with name " << name << std::endl;
    exit(1);
  }
}

void vtkMesh::getPointDataArray(int arrayId, std::vector<double>& data)
{
  if (arrayId < dataSet->GetPointData()->GetNumberOfArrays())
  {
    vtkSmartPointer<vtkDataArray> pd
      = dataSet->GetPointData()->GetArray(arrayId);
    if (pd->GetNumberOfComponents() > 1)
    {
      std::cerr << __func__ << " is only suitable for scalar data, i.e. 1 component\n";
      exit(1);
    }
    data.resize(pd->GetNumberOfTuples());
    double x[1];
    for (int i = 0; i < pd->GetNumberOfTuples(); ++i)
    {
      pd->GetTuple(i,x);
      data[i] = x[0];
    }
  }
  else
  {
    std::cerr << "arrayId exceeds number of point data arrays " << std::endl;
    exit(1);
  }
}

int vtkMesh::getCellDataIdx(const char* name)
{
  // check physical group exist and obtain id
  int idx; 
  dataSet->GetCellData()->GetArray(&name[0u],idx);
  return idx;
}

void vtkMesh::getCellDataArray(const std::string& name, std::vector<double>& data)
{
  int idx;  
  vtkSmartPointer<vtkDataArray> cd 
    = dataSet->GetCellData()->GetArray(&name[0u],idx);
  if (idx != -1)
  {
    if (cd->GetNumberOfComponents() > 1)
    {
      std::cerr << __func__ << " is only suitable for scalar data, i.e. 1 component\n";
      exit(1);
    }
    data.resize(cd->GetNumberOfTuples());
    double x[1];
    for (int i = 0; i < cd->GetNumberOfTuples(); ++i)
    {
      cd->GetTuple(i,x);
      data[i] = x[0];
    }
  }
  else
  {
    std::cerr << "could not find data with name " << name << std::endl;
    exit(1);
  }
}

void vtkMesh::getCellDataArray(int arrayId, std::vector<double>& data)
{
  if (arrayId < dataSet->GetCellData()->GetNumberOfArrays())
  {
    vtkSmartPointer<vtkDataArray> cd
      = dataSet->GetCellData()->GetArray(arrayId);
    if (cd->GetNumberOfComponents() > 1)
    {
      std::cerr << __func__ << " is only suitable for scalar data, i.e. 1 component\n";
      exit(1);
    }
    data.resize(cd->GetNumberOfTuples());
    double x[1];
    for (int i = 0; i < cd->GetNumberOfTuples(); ++i)
    {
      cd->GetTuple(i,x);
      data[i] = x[0];
    }
  }
  else
  {
    std::cerr << "arrayId exceeds number of cell data arrays " << std::endl;
    exit(1);
  }
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



//void addLegacyVTKData(vtkDataArray* arr, const std::string& type, const bool pointOrCell, 
//                      vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp)
//{
//  // casting to appropriate type array
//  if (!type.compare("unsigned_int"))
//  {
//    vtkUnsignedIntArray* typedArr 
//      = vtkUnsignedIntArray::SafeDownCast(arr);
//    pointOrCell ? dataSet_tmp->GetPointData()->AddArray(typedArr)
//                : dataSet_tmp->GetCellData()->AddArray(typedArr);
//  }
//  else if (!type.compare("int"))
//  {
//    std::cout << "in type arr" << std::endl;
//    vtkIntArray* typedArr 
//      = vtkIntArray::SafeDownCast(arr);
//    pointOrCell ? dataSet_tmp->GetPointData()->AddArray(typedArr)
//                : dataSet_tmp->GetCellData()->AddArray(typedArr);
//  }
//  else if (!type.compare("unsigned_long"))
//  {
//    vtkUnsignedLongArray* typedArr 
//      = vtkUnsignedLongArray::SafeDownCast(arr);
//    pointOrCell ? dataSet_tmp->GetPointData()->AddArray(typedArr)
//                : dataSet_tmp->GetCellData()->AddArray(typedArr);
//  }
//  else if (!type.compare("long"))
//  {
//    vtkLongArray* typedArr 
//      = vtkLongArray::SafeDownCast(arr);
//    pointOrCell ? dataSet_tmp->GetPointData()->AddArray(typedArr)
//                : dataSet_tmp->GetCellData()->AddArray(typedArr);
//  }
//  else if (!type.compare("float"))
//  {
//    vtkFloatArray* typedArr 
//      = vtkFloatArray::SafeDownCast(arr);
//    pointOrCell ? dataSet_tmp->GetPointData()->AddArray(typedArr)
//                : dataSet_tmp->GetCellData()->AddArray(typedArr);
//  }
//  else if (!type.compare("double"))
//  {
//    vtkDoubleArray* typedArr 
//      = vtkDoubleArray::SafeDownCast(arr);
//    pointOrCell ? dataSet_tmp->GetPointData()->AddArray(typedArr)
//                : dataSet_tmp->GetCellData()->AddArray(typedArr);
//  }
//}
