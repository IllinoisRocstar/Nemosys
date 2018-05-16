#include <meshBase.H>
#include <vtkMesh.H>
#include <meshGen.H>
#include <TransferBase.H>
#include <SizeFieldBase.H>
#include <Refine.H>
#include <MeshQuality.H>
#include <Cubature.H>
// netgen
namespace nglib
{
  #include<nglib.h>
}
// simmetrix
#ifdef HAVE_SYMMX
  #include <symmxGen.H>
#endif

// TODO: Stop using setPoint/CellDataArray in export methods
//        - instead, use the faster vtkDataArray creation and insertion

meshBase* meshBase::Create(std::string fname)
{
  if (fname.find(".vt") != -1 || fname.find(".stl") != -1 ) 
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

  else if (fname.find(".pntmesh") != -1)
  {
    std::cout << "Detected file in PNTmesh format" << std::endl;
    std::cout << "Processing the file ...." << std::endl;
    return exportPntToVtk(fname);
  }

  else
  {
    int dot = fname.find_last_of('.');
    std::cout << "mesh files with extension " 
              << fname.substr(dot) << " are not supported!" << std::endl;
    exit(1);
  }
  
}

meshBase* meshBase::Create(vtkSmartPointer<vtkDataSet> other, std::string newname)
{
  if (other)
  {
    vtkMesh* vtkmesh = new vtkMesh();
    vtkmesh->dataSet = other;
    vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
    vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();
    vtkmesh->setFileName(newname); 
    return vtkmesh;
  }
  else
  {
    std::cerr << "Nothing to copy!" << std::endl;
    exit(1);
  }
}

std::shared_ptr<meshBase> meshBase::CreateShared(std::string fname)
{
  std::shared_ptr<meshBase> mesh;
  mesh.reset(meshBase::Create(fname));
  return mesh;
}

std::unique_ptr<meshBase> meshBase::CreateUnique(std::string fname)
{
  return std::unique_ptr<meshBase>(meshBase::Create(fname));
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
    generator=nullptr;
  }
  if(!status)
  {
    if (meshEngine == "netgen") 
    {
      std::string newname = trim_fname(fname,".vol");
      return exportVolToVtk(newname);    
    }
    else if (meshEngine == "simmetrix")
    {
      std::string newname = trim_fname(fname, ".vtu");
      return Create(generator->getDataSet(), newname); 
    }
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
  std::unique_ptr<TransferBase> transobj = TransferBase::CreateUnique(method,this,target);
  transobj->setCheckQual(checkQuality);
  return transobj->transferPointData(arrayIDs,newArrayNames);
}

int meshBase::transfer(meshBase* target, std::string method, 
                       const std::vector<std::string>& arrayNames)
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
  std::unique_ptr<TransferBase> transobj = TransferBase::CreateUnique(method,this,target);
  transobj->setCheckQual(checkQuality);
  return transobj->run(newArrayNames); 
}

std::vector<std::vector<double>> 
meshBase::integrateOverMesh(const std::vector<int>& arrayIDs)
{
  std::unique_ptr<GaussCubature> cubature = GaussCubature::CreateUnique(this, arrayIDs);
  return cubature->integrateOverAllCells(); 
}

void meshBase::generateSizeField(std::string method, int arrayID, double dev_mult, bool maxIsmin, double sizeFactor)
{
  std::cout << "Size Factor = " << sizeFactor << std::endl;
  std::unique_ptr<SizeFieldBase> sfobj = SizeFieldBase::CreateUnique(this,method,arrayID,dev_mult,maxIsmin,sizeFactor);
  sfobj->setSizeFactor(sizeFactor);
  sfobj->computeSizeField(arrayID);
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
  //vtkmesh->dataSet = dataSet_tmp->NewInstance();
  //vtkmesh->dataSet->DeepCopy(dataSet_tmp);
  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();
  
  for (int i = 0; i < pointData.size(); ++i)
    vtkmesh->setPointDataArray(&(pointDataNames[i])[0u], pointData[i]);
  for (int i = 0; i < cellData.size(); ++i)
    vtkmesh->setCellDataArray(&(cellDataNames[i])[0u], cellData[i]); 
 
  // Temporary changed to legacy vtk for AMR demos
  // switch back to vtu when done
  std::cout << __FILE__ << __LINE__ << std::endl;
  vtkmesh->setFileName(trim_fname(fname,".vtk"));
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
  //vtkmesh->dataSet = dataSet_tmp->NewInstance();//
  //vtkmesh->dataSet->DeepCopy(dataSet_tmp);//vtkDataSet::SafeDownCast(dataSet_tmp));
  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();

  vtkmesh->setFileName(trim_fname(fname, ".vtu"));
  std::cout << "vtkMesh constructed" << std::endl;

  if(Ngmesh) nglib::Ng_DeleteMesh(Ngmesh);
  nglib::Ng_Exit();
  return vtkmesh;
  
}

// exports pntMesh to vtk format
meshBase* meshBase::exportPntToVtk(std::string fname)
{
  PNTMesh::pntMesh* pMesh;
  pMesh = new PNTMesh::pntMesh(fname);

  vtkMesh* vtkmesh = new vtkMesh();
  
  /*
  if (!pMesh->isCompatible())
  {
    std::cerr << "Mesh contains unsupported element types.\n";
    std::cerr << "Only 3-Node TRIs and 4-Node TETs are supported.\n";
    vtkmesh->numCells = 0;
    vtkmesh->numPoints = 0;
    return vtkmesh;
  }
  */

  // declare points to be pushed into dataSet_tmp
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  // declare dataSet_tmp which will be associated to output vtkMesh
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = 
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  
  int numPoints = pMesh->getNumberOfPoints();
  int numVolCells =  pMesh->getNumberOfCells();

  // allocate memory for points
  points->SetNumberOfPoints(numPoints);
  for (int i = 0; i < numPoints; ++i)
  {
    std::vector<double> point;
    point = pMesh->getPointCrd(i);
    // insert point i
    points->SetPoint(i,&point[0]);     
  }
  // inserting point array into dataSet_tmp
  dataSet_tmp->SetPoints(points);

  // allocating space for cell connectivities
  dataSet_tmp->Allocate(numVolCells);
  for (int i = 0; i < numVolCells; ++i)
  {
    std::cout << "Element " << i << std::endl;
    vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
    std::vector<int> conn;
    VTKCellType vct= 
        pMesh->getVtkCellTag(pMesh->getElmType(i), pMesh->getElmOrder(i));
    conn = pMesh->getElmConn(i, vct);
    for (int j = 0; j < conn.size(); ++j)
    {
      // insert connectivies for cell into cellIds container
      vtkcellIds->InsertNextId(conn[j]);
    }
    // insert connectivies
    dataSet_tmp->InsertNextCell(vct, vtkcellIds); 
  }
  
  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();

  std::cout << "Trimmed name = "
      << trim_fname(fname, ".vtu") << std::endl;
  vtkmesh->setFileName(trim_fname(fname, ".vtu"));
  //vtkmesh->write();
  std::cout << "vtkMesh constructed" << std::endl;

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
                          double edge_scale, std::string ofname, bool transferData,
                          double sizeFactor)
{
  std::unique_ptr<Refine> refineobj
    = std::unique_ptr<Refine>(new Refine(this,method,arrayID,dev_mult,maxIsmin,edge_scale,ofname,sizeFactor));
  refineobj->run(transferData);
}

void meshBase::refineMesh(std::string method, int arrayID, int _order, 
                          std::string ofname, bool transferData)
{
  setOrder(_order); 
  std::unique_ptr<Refine> refineobj
    = std::unique_ptr<Refine>(new Refine(this,method,arrayID,0,0,0,ofname));
  refineobj->run(transferData);
}

void meshBase::refineMesh(std::string method, std::string arrayName, 
                          double dev_mult, bool maxIsmin, 
                          double edge_scale, std::string ofname, bool transferData,
                          double sizeFactor)
{
  int arrayID = IsArrayName(arrayName);
  if (arrayID == -1)
  {
    std::cout << "Array " << arrayName
              << " not fuond in set of point data arrays" << std::endl;
    exit(1);
  }
  refineMesh(method, arrayID, dev_mult, maxIsmin, edge_scale, ofname, transferData, sizeFactor);
}

void meshBase::refineMesh(std::string method, std::string arrayName, int order, 
                          std::string ofname, bool transferData)
{
  int arrayID = IsArrayName(arrayName);
  if (arrayID == -1)
  {
    std::cout << "Array " << arrayName
              << " not fuond in set of point data arrays" << std::endl;
    exit(1);
  }

  refineMesh(method, arrayID, order, ofname, transferData);
}

// added for uniform refinement by driver
void meshBase::refineMesh(std::string method, double edge_scale, std::string ofname, bool transferData)
{
  refineMesh(method, 0, 0, 0, edge_scale, ofname, transferData);
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
  std::unique_ptr<MeshQuality> qualCheck
    = std::unique_ptr<MeshQuality>(new MeshQuality(this)); 
  qualCheck->checkMesh(ofname);
}


int diffMesh(meshBase* mesh1, meshBase* mesh2)
{
  double tol = 1e-6;

  if (mesh1->getNumberOfPoints() != mesh2->getNumberOfPoints() ||
      mesh1->getNumberOfCells() != mesh2->getNumberOfCells())
  {
    std::cerr << "Meshes don't have the same number of points or cells" << std::endl;
    return 1;
  }

  std::cout << mesh1->getNumberOfPoints() << std::endl;
  for (int i = 0; i < mesh1->getNumberOfPoints(); ++i)
  {
    std::vector<double> coord1 = mesh1->getPoint(i);
    std::vector<double> coord2 = mesh2->getPoint(i);
    for (int j = 0; j < 3; ++j)
    {
      std::cout << coord1[j] << std::endl;
      if (std::fabs(coord1[j]-coord2[j]) > tol)
      {
        std::cerr << "Meshes differ in point coordinates" << std::endl;
        return 1;
      }
    } 
  }

  for (int i = 0; i < mesh1->getNumberOfCells(); ++i)
  {
    std::vector<std::vector<double>> cell1 = mesh1->getCellVec(i);
    std::vector<std::vector<double>> cell2 = mesh2->getCellVec(i);
    if (cell1.size() != cell2.size())
    {
      std::cerr << "Meshes differ in cells" << std::endl;
      return 1; 
    }
    for (int j = 0; j < cell1.size(); ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        if (std::fabs(cell1[j][k] - cell2[j][k]) > tol)
        {
          std::cerr << "Meshes differ in cells" << std::endl;
          return 1; 
        }
      }
    }   
  }

  vtkSmartPointer<vtkPointData> pd1 = vtkSmartPointer<vtkPointData>::New();
  pd1 = mesh1->getDataSet()->GetPointData();
  vtkSmartPointer<vtkPointData> pd2 = vtkSmartPointer<vtkPointData>::New();
  pd2 = mesh2->getDataSet()->GetPointData(); 
  int numArr1 = pd1->GetNumberOfArrays(); 
  int numArr2 = pd2->GetNumberOfArrays(); 
 
  if (numArr1 != numArr2)
  {
    std::cerr << "Meshes have different numbers of point data" << std::endl;
    return 1;
  }

  for (int i = 0; i < numArr1; ++i)
  {
    vtkDataArray* da1 = pd1->GetArray(i);
    vtkDataArray* da2 = pd2->GetArray(i);
    int numComponent = da1->GetNumberOfComponents();
    for (int j = 0; j < mesh1->getNumberOfPoints(); ++j)
    {
      double comps1[numComponent];
      double comps2[numComponent];
      da1->GetTuple(j, comps1);
      da2->GetTuple(j, comps2);
      for (int k = 0; k < numComponent; ++k)
      {
        if (std::fabs(comps1[k] - comps2[k]) > tol)
        {
          std::cerr << "Meshes differ in point data values at point " << j <<  std::endl;
          std::cerr << comps1[k] << " " << comps2[k] << std::endl;
          return 1;
        }
      }
    }
  }
  std::cerr << "Meshes are the same" << std::endl;
  return 0;
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



