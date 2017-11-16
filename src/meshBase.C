#include<meshBase.H>
#include<vtkMesh.H>


meshBase* meshBase::Create(std::string fname)
{
  if (fname.find(".vt") != -1 ) 
  {
    vtkMesh* vtkmesh = new vtkMesh(&fname[0u]);
    return vtkmesh;
  }
  
  else if (fname.find(".msh") != -1)
  {
    // convert msh to vtk, trim fname and load
    exportGmshToVtk(fname);
    vtkMesh* vtkmesh = new vtkMesh(&(trim_fname(fname,".vtk"))[0u]);
    return vtkmesh;
  }

  else if (fname.find(".vol") != -1)
  {
    // convert vol to vtk, trim fname and load
  }

  else if (fname.find(".stl") != -1)
  {
    // convert stl to vol to vtk, trim fname and load
  }

  else
  {
    int dot = fname.find_last_of('.');
    std::cout << "mesh files with extension " 
              << fname.substr(dot) << " are not supported!" << std::endl;
    exit(1);
  }
  
}


// get number of points in mesh
int meshUser::getNumberOfPoints() 
{ 
  return mesh->getNumberOfPoints(); 
}

// get number of cells in mesh
int meshUser::getNumberOfCells() 
{ 
  return mesh->getNumberOfCells(); 
}

// get point with id
std::vector<double> meshUser::getPoint(int id) 
{ 
  return mesh->getPoint(id);
}

// get cell with id : returns point indices and respective coordinates
std::map<int,std::vector<double>> meshUser::getCell(int id) 
{ 
  return mesh->getCell(id);
}

// print coordinates of point with id
void meshUser::printPoint(int id)
{ 
  //std::vector<double> coord = mesh->getPoint(id);
  printVec(mesh->getPoint(id));
}

// print point ids and coordinates of cell with id
void meshUser::printCell(int id)
{
  std::map<int, std::vector<double>> cell = mesh->getCell(id);        
  std::map<int, std::vector<double>>::iterator it = cell.begin();
  while(it != cell.end())
  {
    std::cout << it->first << " ";
    printVec(it->second);
    it++;
  }
}


std::pair<int,int> meshBase::exportGmshToVtk(std::string fname)
{
  std::ifstream meshStream(fname);
  if (!meshStream.good())
  {
    std::cout << "Error opening file " << fname << std::endl;
    exit(1);
  }


  std::ofstream vtk(trim_fname(fname,".vtk"));

  if (!vtk.good())
  {
    std::cout << "Error opening: " << fname << std::endl;
    exit(1);
  }

  vtk << "# vtk DataFile Version 2.0" << std::endl 
      << "Converted From Netgen" << std::endl
      << "ASCII" << std::endl
      << "DATASET UNSTRUCTURED_GRID" << std::endl; 

  std::string line;
  int numPoints,numCells, numTri, numTet; 
  numTri = numTet = 0;
  std::vector<std::vector<int>> cells;
  while (getline(meshStream, line))
  {
    if (line.find("$Nodes") != -1)
    {
      getline(meshStream,line);
      std::stringstream ss(line); 
      ss >> numPoints;
      vtk << "POINTS " << numPoints << " double" << std::endl;
      int id;
      double x,y,z;
      for (int i = 0; i < numPoints; ++i)
      {
        getline(meshStream,line);
        std::stringstream ss(line);
        ss >> id >> x >> y >> z;
        vtk << x << " " << y << " " << z << std::endl;          
      }
    }
 

    if (line.find("$Elements") != -1)
    {
      getline(meshStream,line);
      std::stringstream ss(line);
      ss >> numCells;
      int id, type, numTags;
      for (int i = 0; i < numCells; ++i)
      {
        getline(meshStream, line);
        std::stringstream ss(line);
        ss >> id >> type >> numTags;
        std::vector<int> cellIds;
        if (type == 2)
        {
          numTri+=1;
          int tmp;
          for (int j = 0; j < numTags; ++j)
            ss >> tmp;
          for (int j = 0; j < 3; ++j)
          {
            ss >> tmp;
            cellIds.push_back(tmp);
          } 
          cells.push_back(cellIds);
        } 
        else if (type == 4)
        {
          numTet+=1;
          int tmp;
          for (int j = 0; j < numTags; ++j)
            ss >> tmp;
          for (int j = 0; j < 4; ++j)
          {
            ss >> tmp;
            cellIds.push_back(tmp);
          } 
          cells.push_back(cellIds);
        }
        else
        {
          std::cout << "Only triangular and tetrahedral elements are supported!" << std::endl;
          exit(1);
        }
      }
      break;
    }
  }  
  
  vtk << "\nCELLS " << cells.size() << " " << numTri*4 + numTet*5 << std::endl;
  for (int i = 0; i < cells.size(); ++i)
  {
    int size = cells[i].size();
    if (size == 3)
      vtk << 3 << " ";
    else if (size == 4)
      vtk << 4 << " ";
    for (int j = 0; j < size; ++j)
      vtk << cells[i][j] - 1 << " ";
    vtk << std::endl; 
  }
  
  
  vtk << "\nCELL_TYPES " << cells.size() << std::endl;
  for (int i = 0; i < cells.size(); ++i)
  { 
    int size = cells[i].size();
    if (size == 3)
      vtk << VTK_TRIANGLE << std::endl; 
    else if (size == 4)
      vtk << VTK_TETRA << std::endl;
  }

  //TODO: ADD POINT AND CELL DATA AS WELL 
 
  return std::pair<int,int> (numPoints,numTri+numTet);

}



// --------------- AUXILIARY FUNCTIONS ----------------//
template<typename T>
void printVec(const std::vector<T>& v)
{
  for (int i = 0; i < v.size(); ++i)
    std::cout << v[i] << " ";
  std::cout << std::endl;
}

// string trimming for consistent file names 
/*std::string trim_fname(std::string fname, std::string ext)
{
  size_t beg = 0;
  size_t end = fname.find(".");
  std::string name;
  if (end != -1)
  {
    name = fname.substr(beg,end);
    name.append(ext);
    return name;
  }

  else 
  {
    std::cout << "Error finding file extension for " << fname << std::endl;
    exit(1);
  }
} */ 
