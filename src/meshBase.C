#include <meshBase.H>
#include <vtkMesh.H>
#include <meshGen.H>
#include <TransferBase.H>
#include <SizeFieldBase.H>
#include <Refine.H>
#include <MeshQuality.H>

meshBase* meshBase::Create(std::string fname)
{
  if (fname.find(".vt") != -1 ) 
  {
    vtkMesh* vtkmesh = new vtkMesh(&fname[0u]);
    vtkmesh->setFileName(fname);
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

meshBase* meshBase::generateMesh(std::string fname, std::string meshEngine,
                                 meshingParams* params)
{

  if (fname.find(".stl") == -1)
  {
    std::cout << "Only CAD files in STL format are supported" << std::endl;
    exit(1);
  }

  meshGen* generator = meshGen::Create(fname, meshEngine, params);
  int status = generator->createMeshFromSTL(&fname[0u]);
  if (generator) 
  { 
    delete generator;
    generator=0;
  }
  if(!status) 
  {
    std::string newname = trim_fname(fname,".vol");
    return exportVolToVtk(newname);    
    
  }
}


// check for named array in vtk 
int meshBase::IsArrayName(std::string name)
{
  vtkPointData* pd = dataSet->GetPointData();
  if (pd->GetNumberOfArrays()) {
    for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
      std::string curr_name = (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL");
      if (!name.compare(curr_name)) {
        return i;
      }
    }
  }
  return -1;
}


// transfer point data with given ids from this mesh to target
int meshBase::transfer(meshBase* target, std::string method, 
                       const std::vector<int>& arrayIDs)
{
  TransferBase* transobj = TransferBase::Create(method, this, target);//new Transfer(this, target);
  transobj->setCheckQual(checkQuality);
  transobj->runPD(arrayIDs);
  /*for (int i = 0; i < arrayIDs.size(); ++i)
  {
    std::cout << "transferring array " << arrayIDs[i] << std::endl;
    transobj->runPD(arrayIDs[i]); // specify params to run function
  }*/
  if (transobj)
  { 
    delete transobj;
    transobj = 0;
  }
  return 0;
}

int meshBase::transfer(meshBase* target, std::string method, const std::vector<std::string>& arrayNames)
{
  std::vector<int> arrayIDs(arrayNames.size());
  for (int i = 0; i < arrayNames.size(); ++i)
  {
    int id = IsArrayName(arrayNames[i]);
    if (id == -1)
    {
      std::cout << "Array " << arrayNames[i] 
                << " not found in set of point data arrays" << std::endl;
      exit(1);
    }
    arrayIDs[i] = id;
  }
  return transfer(target, method,arrayIDs);

}

// transfer all data from this mesh to target
int meshBase::transfer(meshBase* target, std::string method)
{
  
  TransferBase* transobj = TransferBase::Create(method, this, target);//new Transfer(this, target);
  transobj->setCheckQual(checkQuality);
  int result = transobj->run(); // specify params to run function
  if (transobj)
  { 
    delete transobj;
    transobj = 0;
  }
  return result;
}



void meshBase::generateSizeField(std::string method, int arrayID, double dev_mult, bool maxIsmin)
{
  SizeFieldBase* sfobj = SizeFieldBase::Create(this, method, arrayID, dev_mult, maxIsmin); 
  sfobj->computeSizeField(arrayID);
  if (sfobj)
  {
    delete sfobj;
    sfobj = 0;
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

  std::string line;
  int numPoints,numCells; 
  std::vector<std::vector<std::vector<double>>> pointData;
  std::vector<std::vector<std::vector<double>>> cellData;
  std::vector<std::string> pointDataNames;
  std::vector<std::string> cellDataNames;

  // declare points to be pushed into dataSet_tmp
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  // declare dataSet_tmp which will be associated to output vtkMesh
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = vtkSmartPointer<vtkUnstructuredGrid>::New();
  // map to hold true index of points (gmsh allows non-contiguous ordering)
  std::map<int,int> trueIndex;
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
        trueIndex.insert(std::pair<int,int> (id,i));
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
        vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
        if (type == 2)
        {
          int tmp;
          for (int j = 0; j < numTags; ++j)
            ss >> tmp;
          for (int j = 0; j < 3; ++j)
          {
            ss >> tmp;
            // insert connectivies for cell into cellIds container
            vtkcellIds->InsertNextId(trueIndex[tmp]);//-1);
          }
          // insert connectivies for triangle elements into dataSet 
          dataSet_tmp->InsertNextCell(VTK_TRIANGLE,vtkcellIds); 
        } 
        else if (type == 4)
        {
          int tmp;
          for (int j = 0; j < numTags; ++j)
            ss >> tmp;
          for (int j = 0; j < 4; ++j)
          {
            ss >> tmp;
            // insert connectivities for cell into cellids container
            vtkcellIds->InsertNextId(trueIndex[tmp]);//-1);
          } 
          // insert connectivities for tet elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_TETRA,vtkcellIds);
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
 
  vtkmesh->setFileName(trim_fname(fname,".vtu"));
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
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = 
    vtkSmartPointer<vtkUnstructuredGrid>::New();
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

  vtkmesh->setFileName(trim_fname(fname, ".vtu"));
  std::cout << "vtkMesh constructed" << std::endl;

  if(Ngmesh) nglib::Ng_DeleteMesh(Ngmesh);
  nglib::Ng_Exit();
  return vtkmesh;
  
}


// convert to gmsh format without data
void meshBase::writeMSH(std::ofstream& outputStream)
{
  if(!outputStream.good()) 
  {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
 
  if (!dataSet) 
  {
    std::cout << "No data to write" << std::endl;
    exit(1);
  } 
  // ---------  writing gmsh header ----------- //
  outputStream << "$MeshFormat" << std::endl
               << "2.2 0 8" << std::endl
               << "$EndMeshFormat" << std::endl; 

  // -------- ensure all cell types are tri/tet or below -------------- //
  for (int i = 0; i < numCells; i++)
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
  outputStream << "$Nodes" << std::endl << numPoints << std::endl;
  for (int i = 0; i < numPoints; ++i)
  {
    std::vector<double> pntcrds = getPoint(i);
    outputStream << i + 1 << " "
                 << pntcrds[0] << " "
                 << pntcrds[1] << " "
                 << pntcrds[2] << " " << std::endl;
  }
  outputStream << "$EndNodes" << std::endl;

  // ------------- write element type and connectivity --------------------- //
  outputStream << "$Elements" << std::endl << numCells << std::endl;
  for (int i = 0; i < numCells; ++i)
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
}
// convert to gmsh format with specified point or cell data
void meshBase::writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID)
{

  // write points and cells
  writeMSH(outputStream);
  
  if(!outputStream.good()) 
  {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }
 
  if (!pointOrCell.compare("point"))
  {
    // ---------------------------- write point data ------------------------- //
    vtkPointData* pointData = dataSet->GetPointData(); 
    if (pointData)
    {
      int numArr = pointData->GetNumberOfArrays();
      if (arrayID >= numArr)
      {
        std::cout << "ERROR: arrayID is out of bounds" << std::endl;
        std::cout << "There are " << numArr << " point data arrays" << std::endl;
        exit(1);
      }
      else if (numArr < 1)
      {
        std::cout << "no point data found" << std::endl;
        exit(1);
      }
      vtkDataArray* da = pointData->GetArray(arrayID);
      if(da)
      {
        int numComponent = da->GetNumberOfComponents();
        int numTuple = da->GetNumberOfTuples();
        std::string tmpname = "PointArray";
        tmpname += std::to_string(arrayID);
        outputStream << "$NodeData" << std::endl
                     << 1 << std::endl // 1 string tag
                     << "\"" << (pointData->GetArrayName(arrayID) ? 
                                 pointData->GetArrayName(arrayID) : tmpname) // name of view
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
  else if (!pointOrCell.compare("cell"))
  {
    // -------------------------- write cell data ---------------------------- // 

    vtkCellData* cellData = dataSet->GetCellData();
    if (cellData)
    {
      int numArr = cellData->GetNumberOfArrays();
      if (arrayID >= numArr)
      {
        std::cout << "ERROR: arrayID is out of bounds" << std::endl;
        std::cout << "There are " << numArr << " cell data arrays" << std::endl;
        exit(1);
      }
      else if (numArr < 1)
      {
        std::cout << "no cell data found" << std::endl;
        exit(1);
      }
      vtkDataArray* da = cellData->GetArray(arrayID);
      if (da)
      {
        int numComponent = da->GetNumberOfComponents();
        int numTuple = da->GetNumberOfTuples();
        std::string tmpname = "CellArray";
        tmpname += std::to_string(arrayID);
        outputStream << "$ElementData" << std::endl
                     << 1 << std::endl // 1 string tag
                     << "\"" << (cellData->GetArrayName(arrayID) ? 
                                 cellData->GetArrayName(arrayID) : tmpname) // name of view
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

// convert to gmsh format with specified point or cell data for
void meshBase::writeMSH(std::ofstream& outputStream, std::string pointOrCell, int arrayID, 
                        bool onlyVol)
{
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
  for (int i = 0; i < numCells; i++)
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
  outputStream << "$Nodes" << std::endl << numPoints << std::endl;
  for (int i = 0; i < numPoints; ++i)
  {
    std::vector<double> pntcrds = getPoint(i);
    outputStream << i + 1 << " "
                 << pntcrds[0] << " "
                 << pntcrds[1] << " "
                 << pntcrds[2] << " " << std::endl;
  }
  outputStream << "$EndNodes" << std::endl;

  // ------------- write element type and connectivity --------------------- //
  outputStream << "$Elements" << std::endl << numCells-num_bad << std::endl;
  //int k = 0;
  for (int i = 0; i < numCells; ++i)
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
  vtkCellData* cellData = dataSet->GetCellData();
  vtkDataArray* da = cellData->GetArray(arrayID);
  if (da)
  {
    std::string tmpname = "CellArray";
    tmpname += std::to_string(arrayID);
    outputStream << "$ElementData" << std::endl
                 << 1 << std::endl // 1 string tag
                 << "\"" << (cellData->GetArrayName(arrayID) ? 
                             cellData->GetArrayName(arrayID) : tmpname) // name of view
                 << "\"" << std::endl 
                 << 0 << std::endl // 0 real tag
                 << 3 << std::endl // 3 int tags (dt index, dim of field, number of fields)
                 << 0 << std::endl // dt index
                 << 1 << std::endl // dim of field
                 << numCells-num_bad << std::endl; // number of fields
    for (int j = 0; j < numCells; ++j)
    {
      int type_id = dataSet->GetCellType(j);
      if (type_id == 10)
      {
        double* data = da->GetTuple(j);
        outputStream << j + 1 << " ";
        outputStream << data[0] << " ";
        outputStream << std::endl;
      }
    }
    outputStream << "$EndElementData" << std::endl;    
  }
}

void meshBase::writeMSH(std::string fname, std::string pointOrCell, int arrayID,
                        bool onlyVol)
{
  std::ofstream outputStream(fname.c_str());
  writeMSH(outputStream, pointOrCell, arrayID, onlyVol);
}

void meshBase::writeMSH(std::string fname)
{
  std::ofstream outputStream(fname.c_str());
  writeMSH(outputStream);
}

void meshBase::writeMSH(std::string fname, std::string pointOrCell, int arrayID)
{
  std::ofstream outputStream(fname.c_str());
  writeMSH(outputStream, pointOrCell, arrayID);
}



void meshBase::refineMesh(std::string method, int arrayID, 
                          double dev_mult, bool maxIsmin, 
                          double edge_scale, std::string ofname)
{
  Refine* refineobj = new Refine(this, method, arrayID, dev_mult, maxIsmin, edge_scale, ofname);
  refineobj->run();
  if (refineobj)
  { 
    delete refineobj;
    refineobj = 0;
  }
}

void meshBase::refineMesh(std::string method, std::string arrayName, 
               					  double dev_mult, bool maxIsmin, 
													double edge_scale, std::string ofname)
{
  int arrayID = IsArrayName(arrayName);
  if (arrayID == -1)
  {
    std::cout << "Array " << arrayName
              << " not fuond in set of point data arrays" << std::endl;
    exit(1);
  }
  refineMesh(method, arrayID, dev_mult, maxIsmin, edge_scale, ofname);

}
// added for uniform refinement by driver
void meshBase::refineMesh(std::string method, double edge_scale, std::string ofname)
{
  refineMesh(method, 0, 0, 0, edge_scale, ofname);
}

vtkSmartPointer<vtkCellLocator> meshBase::buildLocator()
{
  vtkSmartPointer<vtkCellLocator> cellLocator = 
    vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(dataSet); 
  cellLocator->BuildLocator();
  return cellLocator;
}

void meshBase::checkMesh(std::string ofname)
{
  MeshQuality* qualcheck = new MeshQuality(this);
  qualcheck->checkMesh(ofname);

  if (qualcheck)
  {
    delete qualcheck;
    qualcheck = 0;
  }

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



