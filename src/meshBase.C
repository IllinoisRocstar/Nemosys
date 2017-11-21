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
//    std::cout << "Detected file in STL format" << std::endl;
//    std::cout << "Generating volume mesh with netgen and exporting to VTK ...." << std::endl;
//    return exportStlToVtk(fname);
  }

  else
  {
    int dot = fname.find_last_of('.');
    std::cout << "mesh files with extension " 
              << fname.substr(dot) << " are not supported!" << std::endl;
    exit(1);
  }
  
}

meshBase* meshBase::generateMesh(std::string fname, std::string meshEngine)
{
  if (meshEngine == "netgen")
  {
    meshNetgen* generator = new meshNetgen();
    int status = generator->createMeshFromSTL(&fname[0u]);
    if (generator) delete generator;
    if(!status) 
    {
      std::string newname = trim_fname(fname,".vol");
      return exportVolToVtk(newname);    
      
    }
  }
}


int meshUser::generateMesh(std::string filename, std::string meshEngine)
{
  if (filename.find(".stl") == -1)
  {
    std::cout << "Only CAD files in STL format are supported" << std::endl;
    exit(1);
  }
  mesh = meshBase::generateMesh(filename,meshEngine);
  write_ext.assign(".vtu");
  fname.assign(trim_fname(filename,".vol"));
  std::cout << "user constructed" << std::endl;
  
  return 0;

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

void meshUser::write()
{
  mesh->write(&(trim_fname(fname,write_ext))[0u], write_ext);
}

meshBase* meshBase::exportGmshToVtk(std::string fname)
{
  std::ifstream meshStream(fname);
  if (!meshStream.good())
  {
    std::cout << "Error opening file " << fname << std::endl;
    exit(1);
  }

  std::string line;
  int numPoints,numCells; 
  std::vector<std::vector<int>> cells;
  std::vector<std::vector<std::vector<double>>> pointData;
  std::vector<std::vector<std::vector<double>>> cellData;
  std::vector<std::string> pointDataNames;
  std::vector<std::string> cellDataNames;

  // declare points to be pushed into dataSet_tmp
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  // declare dataSet_tmp which will be associated to output vtkMesh
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = vtkSmartPointer<vtkUnstructuredGrid>::New();
  while (getline(meshStream, line))
  {
    if (line.find("$Nodes") != -1)
    {
      getline(meshStream,line);
      std::stringstream ss(line); 
      ss >> numPoints;
      int id;
      double x,y,z;
      // allocate memory for points
      points->SetNumberOfPoints(numPoints);
      for (int i = 0; i < numPoints; ++i)
      {
        getline(meshStream,line);
        std::stringstream ss(line);
        ss >> id >> x >> y >> z;
        double point[3];
        point[0] = x; point[1] = y; point[2] = z;
        // insert point i
        points->SetPoint(i,point);     
      }
      // inserting point array into dataSet_tmp
      dataSet_tmp->SetPoints(points);

    }
 

    if (line.find("$Elements") != -1)
    {
      getline(meshStream,line);
      std::stringstream ss(line);
      ss >> numCells;
      int id, type, numTags;
      // allocate space for cell connectivities
      dataSet_tmp->Allocate(numCells);
      for (int i = 0; i < numCells; ++i)
      {
        getline(meshStream, line);
        std::stringstream ss(line);
        ss >> id >> type >> numTags;
        std::vector<int> cellIds;
        vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
        if (type == 2)
        {
          int tmp;
          for (int j = 0; j < numTags; ++j)
            ss >> tmp;
          for (int j = 0; j < 3; ++j)
          {
            ss >> tmp;
            cellIds.push_back(tmp);
            // insert connectivies for cell into cellIds container
            vtkcellIds->InsertNextId(tmp-1);
          }
          // insert connectivies for triangle elements into dataSet 
          dataSet_tmp->InsertNextCell(VTK_TRIANGLE,vtkcellIds); 
          cells.push_back(cellIds);
        } 
        else if (type == 4)
        {
          int tmp;
          for (int j = 0; j < numTags; ++j)
            ss >> tmp;
          for (int j = 0; j < 4; ++j)
          {
            ss >> tmp;
            cellIds.push_back(tmp);
            // insert connectivities for cell into cellids container
            vtkcellIds->InsertNextId(tmp-1);
          } 
          // insert connectivities for tet elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_TETRA,vtkcellIds);
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

  vtkMesh* vtkmesh = new vtkMesh();
  vtkmesh->dataSet = dataSet_tmp->NewInstance();//
  vtkmesh->dataSet->DeepCopy(dataSet_tmp);//vtkDataSet::SafeDownCast(dataSet_tmp));
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();
  
  for (int i = 0; i < pointData.size(); ++i)
    vtkmesh->setPointDataArray(&(pointDataNames[i])[0u], pointData[i]);
  for (int i = 0; i < cellData.size(); ++i)
    vtkmesh->setCellDataArray(&(cellDataNames[i])[0u], cellData[i]); 
 

  std::cout << "vtkMesh constructed" << std::endl;

  return vtkmesh;

}

meshBase* meshBase::exportVolToVtk(std::string fname)
{
  nglib::Ng_Mesh* Ngmesh;
  nglib::Ng_Init();
  Ngmesh = nglib::Ng_NewMesh();
  
  int status = nglib::Ng_MergeMesh(Ngmesh, &fname[0u]);
  if (status)
  {
    std::cout << "Error: NetGen Returned: " << status << std::endl;
    std::cout << "Could not load " << fname << " into netgen" << std::endl;
    exit(1); 
  } 

  // declare points to be pushed into dataSet_tmp
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  // declare dataSet_tmp which will be associated to output vtkMesh
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = vtkSmartPointer<vtkUnstructuredGrid>::New();
  int numNgPoints = nglib::Ng_GetNP(Ngmesh);
  int numSurfCells = nglib::Ng_GetNSE(Ngmesh); 
  int numVolCells = nglib::Ng_GetNE(Ngmesh);

  // allocate memory for points
  points->SetNumberOfPoints(numNgPoints);
  for (int i = 1; i <= numNgPoints; ++i)
  {
    double point[3];
    nglib::Ng_GetPoint(Ngmesh,i,point);
    // insert point i
    points->SetPoint(i-1,point);     
  }
  // inserting point array into dataSet_tmp
  dataSet_tmp->SetPoints(points);

  // allocating space for cell connectivities
  dataSet_tmp->Allocate(numSurfCells + numVolCells);

  for (int i = 1; i <= numSurfCells; ++i)
  {
    vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
    int tri[3];
    nglib::Ng_GetSurfaceElement(Ngmesh,i,tri);  
    for (int j = 0; j < 3; ++j)
    {
      // insert connectivies for cell into cellIds container
      vtkcellIds->InsertNextId(tri[j]-1);
    }
    // insert connectivies for triangle elements into dataSet 
    dataSet_tmp->InsertNextCell(VTK_TRIANGLE,vtkcellIds); 
  }

  for (int i = 1; i <= numVolCells; ++i)
  { 
    vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
    int tet[4];
    nglib::Ng_GetVolumeElement(Ngmesh, i, tet);
    for (int j = 0; j < 4; ++j)
    {
      // insert connectivies for cell into cellIds container
      vtkcellIds->InsertNextId(tet[j]-1);
    }
    // insert connectivies for triangle elements into dataSet 
    dataSet_tmp->InsertNextCell(VTK_TETRA,vtkcellIds); 

  }
  
  vtkMesh* vtkmesh = new vtkMesh();
  vtkmesh->dataSet = dataSet_tmp->NewInstance();//
  vtkmesh->dataSet->DeepCopy(dataSet_tmp);//vtkDataSet::SafeDownCast(dataSet_tmp));
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();

  std::cout << "vtkMesh constructed" << std::endl;

  if(Ngmesh) nglib::Ng_DeleteMesh(Ngmesh);
  nglib::Ng_Exit();
  return vtkmesh;
  
}



/*meshBase* meshBase::exportStlToVtk(std::string fname)
{
  netgenInterface* tmp = new netgenInterface();
  tmp->createMeshFromSTL(&fname[0u]);
  char* name = &(trim_fname(fname,".vtk"))[0u];
  tmp->exportToVTK(name);
  delete tmp;
  vtkMesh* vtkmesh = new vtkMesh(name);
  return vtkmesh;
}*/



