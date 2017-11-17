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
    std::cout << "Detected file in GMSH format" << std::endl;
    std::cout << "Exporting to VTK ...." << std::endl;
    return exportGmshToVtk(fname);
  }

  else if (fname.find(".vol") != -1)
  {
    std::cout << "Detected file in Netgen .vol format" << std::endl;
    std::cout << "Exporting to VTK ...." << std::endl;
    return exportVolToVtk(fname);
  }

  else if (fname.find(".stl") != -1)
  {
    std::cout << "Detected file in STL format" << std::endl;
    std::cout << "Generating volume mesh with netgen and exporting to VTK ...." << std::endl;
    return exportStlToVtk(fname);
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

void meshUser::report()
{
  mesh->report(&fname[0u]);
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

meshBase* meshBase::exportGmshToVtk(std::string fname)
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

  vtk << "# vtk DataFile Version 3.0" << std::endl 
      << "Converted From GMSH MSH" << std::endl
      << "ASCII" << std::endl
      << "DATASET UNSTRUCTURED_GRID" << std::endl; 

  std::string line;
  int numPoints,numCells, numTri, numTet; 
  numTri = numTet = 0;
  std::vector<std::vector<int>> cells;
  std::vector<std::vector<std::vector<double>>> pointData;
  std::vector<std::vector<std::vector<double>>> cellData;
  std::vector<std::string> pointDataNames;
  std::vector<std::string> cellDataNames;
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
    }
  
    if (line.find("$NodeData") != -1)
    {
      std::vector<std::vector<double>> currPointData;
      getline(meshStream,line); // moving to num_string_tags
      std::stringstream ss(line);
      int num_string_tags;
      ss >> num_string_tags;
      std::string dataname;
      for (int i = 0; i < num_string_tags; ++i)
      { 
        getline(meshStream,line); // get string tag
        if (i == 0)
        {
          std::stringstream ss(line);
          ss >> dataname;
        }
      }
      dataname.erase(std::remove(dataname.begin(),dataname.end(),'\"'),dataname.end());
      pointDataNames.push_back(dataname);
      getline(meshStream,line); // moving to num_real_tags
      std::stringstream ss1(line);
      int num_real_tags;
      ss1 >> num_real_tags;
      for (int i = 0; i < num_real_tags; ++i)
        getline(meshStream,line);
 
      getline(meshStream,line); // moving to num_int_tags
      std::stringstream ss2(line);
      int num_int_tags;
      ss2 >> num_int_tags;
      int dt,dim,numFields,tmp;
      for (int i = 0; i < num_int_tags; ++i)
      {   
        getline(meshStream,line); // get int tag
        std::stringstream ss(line);
        if (i == 0)
          ss >> dt;
        else if (i == 1)
          ss >> dim;
        else if (i == 2)
          ss >> numFields;
        else
          ss >> tmp;
      }
      for (int i = 0; i < numFields; ++i) 
      { 
        std::vector<double> data(dim);
        int id;
        double val;
        getline(meshStream,line);
        std::stringstream ss(line);
        ss >> id;
        for (int j = 0; j < dim; ++j)
        {
          ss >> val;
          data[j] = val;
        }
        currPointData.push_back(data);
      }
      pointData.push_back(currPointData);
    }
    
    if (line.find("$ElementData") != -1)
    {
      std::vector<std::vector<double>> currCellData;
      getline(meshStream,line); // moving to num_string_tags
      std::stringstream ss(line);
      int num_string_tags;
      ss >> num_string_tags;
      std::string dataname;
      for (int i = 0; i < num_string_tags; ++i)
      { 
        getline(meshStream,line); // get string tag
        if (i == 0)
        {
          std::stringstream ss(line);
          ss >> dataname;
        }
      }
      dataname.erase(std::remove(dataname.begin(),dataname.end(),'\"'),dataname.end());
      cellDataNames.push_back(dataname);
      getline(meshStream,line); // moving to num_real_tags
      std::stringstream ss1(line);
      int num_real_tags;
      ss1 >> num_real_tags;
      for (int i = 0; i < num_real_tags; ++i)
        getline(meshStream,line);
 
      getline(meshStream,line); // moving to num_int_tags
      std::stringstream ss2(line);
      int num_int_tags;
      ss2 >> num_int_tags;
      int dt,dim,numFields,tmp;
      for (int i = 0; i < num_int_tags; ++i)
      {   
        getline(meshStream,line); // get int tag
        std::stringstream ss(line);
        if (i == 0)
          ss >> dt;
        else if (i == 1)
          ss >> dim;
        else if (i == 2)
          ss >> numFields;
        else
          ss >> tmp;
      }
      for (int i = 0; i < numFields; ++i) 
      { 
        std::vector<double> data(dim);
        int id;
        double val;
        getline(meshStream,line);
        std::stringstream ss(line);
        ss >> id;
        for (int j = 0; j < dim; ++j)
        {
          ss >> val;
          data[j] = val;
        }
        currCellData.push_back(data);
      }
      cellData.push_back(currCellData);
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

  vtk.close();
  vtkMesh* vtkmesh = new vtkMesh(&(trim_fname(fname,".vtk"))[0u]);
  for (int i = 0; i < pointData.size(); ++i)
    vtkmesh->setPointDataArray(&(pointDataNames[i])[0u], pointData[i]);
  for (int i = 0; i < cellData.size(); ++i)
    vtkmesh->setCellDataArray(&(cellDataNames[i])[0u], cellData[i]); 

  return vtkmesh;

}

meshBase* meshBase::exportVolToVtk(std::string fname)
{
  netgenInterface* tmp = new netgenInterface();
  int status = tmp->importFromVol(&fname[0u]);
  delete tmp;
  if (!status)
  {
    char* name = &(trim_fname(fname,".vtk"))[0u];
    tmp->exportToVTK(name);
    vtkMesh* vtkmesh = new vtkMesh(name);
    return vtkmesh;
  }
  else
  {
    std::cout << "error converting file to vtk" << std::endl;
    exit(1);
  }
}

meshBase* meshBase::exportStlToVtk(std::string fname)
{
  netgenInterface* tmp = new netgenInterface();
  tmp->createMeshFromSTL(&fname[0u]);
  char* name = &(trim_fname(fname,".vtk"))[0u];
  tmp->exportToVTK(name);
  delete tmp;
  vtkMesh* vtkmesh = new vtkMesh(name);
  return vtkmesh;
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
