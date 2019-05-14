#include <meshBase.H>
#include <vtkMesh.H>
#include <meshGen.H>
#include <TransferBase.H>
#include <SizeFieldBase.H>
#include <Refine.H>
#include <MeshQuality.H>
#include <Cubature.H>
#include <meshPartitioner.H>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkCellTypes.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkAppendFilter.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>

// netgen
namespace nglib
{
  #include<nglib.h>
}
// simmetrix
#ifdef HAVE_SYMMX
  #include <symmxGen.H>
#endif

// stl
#include <algorithm>

// aux
#include <AuxiliaryFunctions.H>

#ifdef HAVE_EXODUSII
#include "exodusII.h"
#include "exoMesh.H"
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

  else if (fname.find(".g") != -1 || fname.find(".exo") != -1)
  {
    std::cout << "Detected file in Exodus II format" << std::endl;
    std::cout << "Processing the file ...." << std::endl;
    return exportExoToVtk(fname);
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
  return new vtkMesh(other, newname);
}

meshBase* meshBase::Create(const std::vector<double>& xCrds,
                           const std::vector<double>& yCrds,
                           const std::vector<double>& zCrds,
                           const std::vector<int>& elmConn, const int cellType,
                           std::string newname)
{
  return new vtkMesh(xCrds, yCrds, zCrds, elmConn, cellType, newname);
}

std::unique_ptr<meshBase>
meshBase::CreateUnique(const std::vector<double>& xCrds,
                       const std::vector<double>& yCrds,
                       const std::vector<double>& zCrds,
                       const std::vector<int>& elmConn, const int cellType,
                       std::string newname)
{
  return std::unique_ptr<meshBase>
          (meshBase::Create(xCrds,yCrds,zCrds,elmConn,cellType,newname));
}

std::unique_ptr<meshBase> 
meshBase::CreateUnique(vtkSmartPointer<vtkDataSet> other, std::string newname)
{
  return std::unique_ptr<meshBase>(meshBase::Create(other,newname));
}

std::unique_ptr<meshBase>
meshBase::CreateUnique(meshBase* mesh)
{
  return std::unique_ptr<meshBase>(mesh);
}

std::shared_ptr<meshBase>
meshBase::CreateShared(meshBase* mesh)
{
  std::shared_ptr<meshBase> sharedMesh;
  sharedMesh.reset(mesh);
  return sharedMesh;
} 

std::shared_ptr<meshBase> 
meshBase::CreateShared(vtkSmartPointer<vtkDataSet> other, std::string newname)
{
  std::shared_ptr<meshBase> mesh;
  mesh.reset(meshBase::Create(other,newname));
  return mesh;
}
std::shared_ptr<meshBase> meshBase::CreateShared(std::string fname)
{
  std::shared_ptr<meshBase> mesh;
  mesh.reset(meshBase::Create(fname));
  return mesh;
}

std::shared_ptr<meshBase> meshBase::CreateShared(const std::vector<double>& xCrds,
                                                 const std::vector<double>& yCrds,
                                                 const std::vector<double>& zCrds,
                                                 const std::vector<int>& elmConn, 
                                                 const int cellType,
                                                 std::string newname)
{
  std::shared_ptr<meshBase> mesh;
  mesh.reset(meshBase::Create(xCrds, yCrds, zCrds, elmConn, cellType, newname));
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
    std::cerr << "Only CAD files in STL format are supported" << std::endl;
    exit(1);
  }

  meshGen* generator = meshGen::Create(fname, meshEngine, params);
  if (generator)
  {
    int status = generator->createMeshFromSTL(&fname[0u]);
    meshBase* ret;
    if(!status)
    {
      if (meshEngine == "netgen") 
      {
        std::string newname = trim_fname(fname,".vol");
        ret = exportVolToVtk(newname);    
      }
      else if (meshEngine == "simmetrix")
      {
        std::string newname = trim_fname(fname, ".vtu");
        ret = Create(generator->getDataSet(), newname); 
      }
    }
    delete generator;
    generator=nullptr;
    return ret;
  }
  else
  {
    std::cerr << "Could not create mesh generator" << std::endl;
    exit(1);
  }
}

meshBase* meshBase::stitchMB(const std::vector<meshBase*>& mbObjs)
{
  if (mbObjs.size())
  {
    vtkSmartPointer<vtkAppendFilter> appender
      = vtkSmartPointer<vtkAppendFilter>::New();
    appender->MergePointsOn();
    for (int i = 0; i < mbObjs.size(); ++i)
    {
      appender->AddInputData(mbObjs[i]->getDataSet());
    }
    appender->Update();
    return meshBase::Create(appender->GetOutput(), "stitched.vtu");
  }
  else
  {
    std::cerr << "Nothing to stitch!" << std::endl;
    exit(1);
  }
}

std::shared_ptr<meshBase> 
meshBase::stitchMB(const std::vector<std::shared_ptr<meshBase>>& _mbObjs)
{
  std::vector<meshBase*> mbObjs(_mbObjs.size());
  for (int i = 0; i < mbObjs.size(); ++i)
  {
    mbObjs[i] = _mbObjs[i].get();
  }
  return meshBase::CreateShared(meshBase::stitchMB(mbObjs));
}

meshBase* meshBase::extractSelectedCells(meshBase* mesh, const std::vector<int>& cellIds)
{
  vtkSmartPointer<vtkIdTypeArray> selectionIds 
    = vtkSmartPointer<vtkIdTypeArray>::New();
  selectionIds->SetNumberOfComponents(1);
  for (int i = 0; i < cellIds.size(); ++i)
  {
    selectionIds->InsertNextValue(cellIds[i]);
  }
  return extractSelectedCells(mesh->getDataSet(), selectionIds);
}

meshBase* meshBase::extractSelectedCells(vtkSmartPointer<vtkDataSet> mesh,
                                         vtkSmartPointer<vtkIdTypeArray> cellIds)
{
  vtkSmartPointer<vtkSelectionNode> selectionNode
    = vtkSmartPointer<vtkSelectionNode>::New(); 
  selectionNode->SetFieldType(vtkSelectionNode::CELL);
  selectionNode->SetContentType(vtkSelectionNode::INDICES);
  selectionNode->SetSelectionList(cellIds);
  // create the selection
  vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
  selection->AddNode(selectionNode);
  // creating extractor
  vtkSmartPointer<vtkExtractSelection> extractSelection
    = vtkSmartPointer<vtkExtractSelection>::New();
  // set stitch surf as input on first port
  extractSelection->SetInputData(0, mesh);
  // set selectionNode as input on second port
  extractSelection->SetInputData(1, selection);
  extractSelection->Update();
  vtkSmartPointer<vtkDataSet> extractedCellMesh
    = vtkDataSet::SafeDownCast(extractSelection->GetOutput());
  meshBase* selectedCellMesh 
    = meshBase::Create(extractedCellMesh, "extracted.vtu");
  return selectedCellMesh;
}

// check for named array in vtk 
int meshBase::IsArrayName(std::string name, const bool pointOrCell)
{
  if (!pointOrCell)
  {
    vtkPointData* pd = dataSet->GetPointData();
    if (pd->GetNumberOfArrays()) 
    {
      for (int i = 0; i < pd->GetNumberOfArrays(); ++i) 
      {
        std::string curr_name = (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL");
        if (!name.compare(curr_name)) 
        {
          return i;
        }
      }
    }
  }
  else 
  {
    vtkCellData* cd = dataSet->GetCellData();
    if (cd->GetNumberOfArrays())
    {
      for (int i = 0; i < cd->GetNumberOfArrays(); ++i)
      {
        std::string curr_name = (cd->GetArrayName(i) ? cd->GetArrayName(i) : "Null");
        if (!name.compare(curr_name))
        {
          return i;
        }
      }
    }
  }
  return -1;
}


// transfer point data or cell data with given ids from this mesh to target
int meshBase::transfer(meshBase* target, std::string method, 
                       const std::vector<int>& arrayIDs, bool pointOrCell)
{
  std::unique_ptr<TransferBase> transobj = TransferBase::CreateUnique(method,this,target);
  transobj->setCheckQual(checkQuality);
  // adding continuity flag, set by this object
  transobj->setContBool(continuous);
  if (!pointOrCell)
  {
    transobj->transferPointData(arrayIDs, newArrayNames);
  }
  else
  {
    transobj->transferCellData(arrayIDs, newArrayNames);
  }
}

int meshBase::transfer(meshBase* target, std::string method, 
                       const std::vector<std::string>& arrayNames, bool pointOrCell)
{
  std::vector<int> arrayIDs(arrayNames.size());
  for (int i = 0; i < arrayNames.size(); ++i)
  {
    int id = IsArrayName(arrayNames[i], pointOrCell);
    if (id == -1)
    {
      std::cout << "Array " << arrayNames[i] 
                << " not found in set of data arrays" << std::endl;
      exit(1);
    }
    arrayIDs[i] = id;
  }
  return transfer(target, method,arrayIDs, pointOrCell);
}

// transfer all data from this mesh to target
int meshBase::transfer(meshBase* target, std::string method)
{
  std::unique_ptr<TransferBase> transobj = TransferBase::CreateUnique(method,this,target);
  transobj->setCheckQual(checkQuality);
  transobj->setContBool(continuous);
  return transobj->run(newArrayNames);
}

// partition mesh into numPartition pieces (static fcn)
std::vector<std::shared_ptr<meshBase>> 
meshBase::partition(const meshBase* mbObj, const int numPartitions)
{
  // construct partitioner with meshBase object
  meshPartitioner* mPart = new meshPartitioner(mbObj);
  if (mPart->partition(numPartitions))
  {
    exit(1); 
  }
  // initialize vector of meshBase partitions
  std::vector<std::shared_ptr<meshBase>> mbParts(numPartitions); 
  for (int i = 0; i < numPartitions; ++i)
  {
    // define coordinates
    std::vector<std::vector<double>> comp_crds(mbObj->getVertCrds()); 
    // get partition connectivity and zero index it
    std::vector<int> vtkConn(mPart->getConns(i));
    for (auto it = vtkConn.begin(); it != vtkConn.end(); ++it)
    {
      *it -= 1;
    }
    std::string basename(trim_fname(mbObj->getFileName(), ""));
    basename += std::to_string(i);
    basename += ".vtu";
    // construct meshBase partition from coordinates and connectivities from partitioner
    mbParts[i] = meshBase::CreateShared(mPart->getCrds(i, comp_crds[0]),
                                        mPart->getCrds(i, comp_crds[1]),
                                        mPart->getCrds(i, comp_crds[2]),
                                        vtkConn, VTK_TETRA, basename);
    // add partition id to node and cell data of mbPart
    vtkSmartPointer<vtkIntArray> nodePartitionIds = vtkSmartPointer<vtkIntArray>::New();
    nodePartitionIds->SetName("NodePartitionIds");
    nodePartitionIds->SetNumberOfComponents(1);
    nodePartitionIds->SetNumberOfTuples(mbParts[i]->getNumberOfPoints());
    nodePartitionIds->FillComponent(0, i);
    mbParts[i]->getDataSet()->GetPointData()->AddArray(nodePartitionIds);

    vtkSmartPointer<vtkIntArray> cellPartitionIds = vtkSmartPointer<vtkIntArray>::New();
    cellPartitionIds->SetName("CellPartitionIds");
    cellPartitionIds->SetNumberOfComponents(1);
    cellPartitionIds->SetNumberOfTuples(mbParts[i]->getNumberOfCells());
    cellPartitionIds->FillComponent(0, i);
    mbParts[i]->getDataSet()->GetCellData()->AddArray(cellPartitionIds);

    // add global node index array to partition
    std::map<int,int> partToGlobNodeMap(mPart->getPartToGlobNodeMap(i)); 
    vtkSmartPointer<vtkIdTypeArray> globalNodeIds = vtkSmartPointer<vtkIdTypeArray>::New();
    globalNodeIds->SetName("GlobalNodeIds");
    globalNodeIds->SetNumberOfComponents(1);
    globalNodeIds->SetNumberOfTuples(mbParts[i]->getNumberOfPoints());
    globalNodeIds->SetNumberOfValues(mbParts[i]->getNumberOfPoints());
    auto it = partToGlobNodeMap.begin();
    int idx = 0;
    while (it != partToGlobNodeMap.end())
    {
      int globidx = it->second-1;
      int locidx = it->first-1;
      globalNodeIds->SetTuple1(idx,globidx);
      mbParts[i]->globToPartNodeMap[globidx] = locidx; 
      mbParts[i]->partToGlobNodeMap[locidx] = globidx;
      ++idx; 
      ++it;
    }
    mbParts[i]->getDataSet()->GetPointData()->AddArray(globalNodeIds);
    // add global cell index array to partition
    std::map<int,int> partToGlobCellMap(mPart->getPartToGlobElmMap(i));
    //vtkSmartPointer<vtkIdTypeArray> globalCellIds = vtkSmartPointer<vtkIdTypeArray>::New();
    //globalCellIds->SetName("GlobalCellIds");
    //globalCellIds->SetNumberOfComponents(1);
    //globalCellIds->SetNumberOfTuples(mbPart->getNumberOfCells());
    //globalCellIds->SetNumberOfValues(mbPart->getNumberOfCells());
    it = partToGlobCellMap.begin();
    idx = 0;
    while (it != partToGlobCellMap.end())
    {
      int globidx = it->second-1;
      int locidx = it->first-1;
  //    globalCellIds->SetTuple1(idx,globidx);
      mbParts[i]->globToPartCellMap[globidx] = locidx;
      mbParts[i]->partToGlobCellMap[locidx] = globidx;
      ++idx;
      ++it;
    }
    //mbPart->getDataSet()->GetCellData()->AddArray(globalCellIds);  
    mbParts[i]->write();
    //mbParts[i] = mbPart;
  }
  delete mPart; mPart = nullptr;
  return mbParts;
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
  std::unique_ptr<SizeFieldBase> sfobj 
    = SizeFieldBase::CreateUnique(this,method,arrayID,dev_mult,maxIsmin,sizeFactor);
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

  bool warning = true;

  std::string line;
  int numPoints,numCells,numPhysGrps;
  bool fndPhyGrp = false;
  std::vector<int> physGrpDims;
  std::map<int, std::string> physGrpIdName;
  std::vector<std::vector<std::vector<double>>> pointData;
  std::vector<std::vector<std::vector<double>>> cellData;
  std::vector<std::vector<double>> cellPhysGrpIds;
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
    if (line.find("$PhysicalNames") != -1)
    {
      fndPhyGrp = true;
      getline(meshStream,line);
      std::stringstream ss(line); 
      ss >> numPhysGrps;
      std::cout << "Found " << numPhysGrps << " physical groups!\n";
      int grpDim, grpId;
      std::string grpName;
      for (int i = 0; i < numPhysGrps; ++i)
      {
        getline(meshStream,line);
        std::stringstream ss(line);
        ss >> grpDim >> grpId >> grpName;
        grpName.erase(std::remove(grpName.begin(),grpName.end(),'\"'),grpName.end());
        //std::cout << "Group Dim = " << grpDim
        //    << "\nGroup Id = " << grpId
        //    << "\nGroup Name = " << grpName 
        //    << std::endl;
        physGrpDims.push_back(grpDim);
        physGrpIdName[grpId] = grpName;
      }
    }
    
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
      //std::cout << "line = " << line << std::endl;
      std::stringstream ss(line);
      ss >> numCells;
      int id, type, numTags;
      double tmp2[1];
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
          if (!fndPhyGrp)
          {
            for (int j = 0; j < numTags; ++j)
                ss >> tmp;
          } 
          else 
          {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags-1; ++j)
                ss >> tmp;
          }
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
          if (!fndPhyGrp)
          {
            for (int j = 0; j < numTags; ++j)
                ss >> tmp;
          } 
          else 
          {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags-1; ++j)
                ss >> tmp;
          }
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
          if (warning)
          {
            std::cout << "Warning: Only triangular and tetrahedral elements are supported, rest is ignored! " << std::endl;
            warning = false;
            //exit(1);
          }
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

  if (fndPhyGrp)
  {
      cellDataNames.push_back("PhysGrpId");
      cellData.push_back(cellPhysGrpIds);
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

// exports exodusII to vtk format
meshBase* meshBase::exportExoToVtk(std::string fname)
{
#ifdef HAVE_EXODUSII
  // opening the file
  int CPU_word_size,IO_word_size;
  int fid, _exErr;
  float version;
  CPU_word_size = sizeof(float);
  IO_word_size = 0;

  /* open EXODUS II files */
  fid = ex_open(fname.c_str(),EX_READ,&CPU_word_size,&IO_word_size,&version);
  EXOMesh::wrnErrMsg(_exErr, "Problem opening file "+fname+"\n");

  // building mesh database
  vtkMesh* vtkmesh = new vtkMesh();

  // declare points to be pushed into dataSet_tmp
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // declare dataSet_tmp which will be associated to output vtkMesh
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp = 
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  
  int numPoints;
  int numVolCells;
  int numElmBlk;
  int numNdeSet;
  int numSideSet;

  // parameter inquiry from Exodus file
  int num_props;
  float fdum;
  char *cdum;
  _exErr = ex_inquire(fid, EX_INQ_API_VERS, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II API version is "<< fdum << std::endl;
  _exErr = ex_inquire(fid, EX_INQ_DB_VERS, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II Database version is "<< fdum << std::endl;
  _exErr = ex_inquire(fid, EX_INQ_DIM, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of coordinate dimensions is "<< num_props << std::endl;
  if (num_props != 3)
    EXOMesh::wrnErrMsg(-1, "Only 3D mesh data is supported!\n");
  _exErr = ex_inquire(fid, EX_INQ_NODES, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numPoints = num_props;
  std::cout << "Number of points "<< numPoints << std::endl;
  _exErr = ex_inquire(fid, EX_INQ_ELEM, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numVolCells =  num_props;
  std::cout << "Number of elements "<< numVolCells << std::endl;
  _exErr = ex_inquire(fid, EX_INQ_ELEM_BLK, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numElmBlk = num_props;
  std::cout << "Number of element blocks "<< numElmBlk << std::endl;
  _exErr = ex_inquire(fid, EX_INQ_NODE_SETS, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numNdeSet = num_props;
  std::cout << "Number of node sets "<< numNdeSet << std::endl;
  _exErr = ex_inquire(fid, EX_INQ_SIDE_SETS, &num_props, &fdum, cdum);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numSideSet = num_props;
  std::cout << "Number of side sets "<< numSideSet << std::endl;

  // nodal coordinates
  std::vector<float> x,y,z;
  x.resize(numPoints,0);
  y.resize(numPoints,0);
  z.resize(numPoints,0);
  _exErr = ex_get_coord(fid, &x[0], &y[0], &z[0]);
  EXOMesh::wrnErrMsg(_exErr, "Problem reading nodal coordinates.\n");
  
  // allocate memory for points
  points->SetNumberOfPoints(numPoints);
  for (int i=0; i<numPoints; ++i)
  {
    std::vector<double> pt = {x[i], y[i], z[i]};
    points->SetPoint(i,&pt[0]); 
  }
  // inserting point array into dataSet_tmp
  dataSet_tmp->SetPoints(points);

  // Connectivities
  // allocating space for cell connectivities
  dataSet_tmp->Allocate(numVolCells);
  // read element blocks
  for (int iEB=1; iEB<=numElmBlk; iEB++)
  {
    int num_el_in_blk, num_nod_per_el, num_attr, *connect;
    float *attrib;
    char elem_type[MAX_STR_LENGTH+1];
    // read element block parameters
    _exErr = ex_get_elem_block (fid, iEB, elem_type, &num_el_in_blk, &num_nod_per_el, &num_attr);
    EXOMesh::wrnErrMsg(_exErr, "Problem reading element block parameters.\n");
    // read element connectivity
    std::vector<int> conn;
    conn.resize(num_el_in_blk*num_nod_per_el,0);
    _exErr = ex_get_elem_conn (fid, iEB, &conn[0]);
    EXOMesh::wrnErrMsg(_exErr, "Problem reading element block connectivites.\n");
    // read element block attributes
    //std::vector<float> attr;
    //attr.resize(0.,num_el_in_blk*num_nod_per_el);
    //_exErr = ex_get_elem_attr (fid, iEB, &attrib[0]);
    //EXOMesh::wrnErrMsg(_exErr, "Problem reading element block attributes.\n");
    for (int iEl = 0; iEl < num_el_in_blk; ++iEl)
    {
      vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
      VTKCellType vct = EXOMesh::e2vEMap(EXOMesh::elmTypeNum(elem_type));
      for (int jc=iEl*num_nod_per_el; 
               jc<(iEl+1)*num_nod_per_el;
               jc++)
      {
        // insert connectivies for cell into cellIds container
        vtkcellIds->InsertNextId(conn[jc]-1);
      }
      // insert connectivies
      dataSet_tmp->InsertNextCell(vct, vtkcellIds); 
    }
  }

  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();

  std::cout << "Trimmed name = "
      << trim_fname(fname, ".vtu") << std::endl;
  vtkmesh->setFileName(trim_fname(fname, ".vtu"));
  //vtkmesh->write();
  std::cout << "vtkMesh constructed" << std::endl;

  // closing the file
  _exErr = ex_close(fid);
  EXOMesh::wrnErrMsg(_exErr, "Problem closing the exodusII file.");
  
  return vtkmesh;
#else
  std::cerr << "Error: Compile with Exodus II to use this.\n";
  exit(-1);
#endif
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

void writePatchMap(const std::string& mapFile, const std::map<int,int>& patchMap)
{
  std::ofstream outputStream(mapFile);
  if (!outputStream.good())
  {
    std::cout << "Error opening file " << mapFile << std::endl;
    exit(1);
  }
  writePatchMap(outputStream, patchMap);
}

void writePatchMap(std::ofstream& outputStream, const std::map<int,int>& patchMap)
{
  outputStream << patchMap.size() << std::endl;
  outputStream << patchMap.size() << std::endl;
  auto it = patchMap.begin();
  int normPatchNo = 1;
  while (it != patchMap.end())
  {
    for (int i = 0; i < 2; ++i)
    {
      outputStream << std::setw(2) << std::left << it->first << " ";
    }
    outputStream << std::setw(2) << std::left << normPatchNo << " ";
    outputStream << std::endl;
    ++it;
    normPatchNo++;
  }
}

void meshBase::writeCobalt(meshBase* surfWithPatches,
                           const std::string& mapFile, std::ofstream& outputStream)
{
  if (!surfWithPatches) 
  {
    std::cout << "surface mesh is empty!" << std::endl;
    exit(1);
  }
  if (surfWithPatches->IsArrayName("patchNo", 1) == -1)
  {
    std::cout << "surface mesh must have patchNo cell array" << std::endl;
    exit(1);
  } 
  vtkSmartPointer<vtkIdList> facePtIds;
  vtkSmartPointer<vtkIdList> sharedCellPtIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 = vtkSmartPointer<vtkGenericCell>::New(); 
  vtkSmartPointer<vtkGenericCell> genCell2 = vtkSmartPointer<vtkGenericCell>::New(); 
  std::map<std::vector<int>, std::pair<int,int>, sortIntVec_compare> faceMap;
  // building cell locator for looking up patch number in remeshed surface mesh
  vtkSmartPointer<vtkCellLocator> surfCellLocator = surfWithPatches->buildLocator(); 
  // maximum number of vertices per face (to be found in proceeding loop)
  int nVerticesPerFaceMax = 0;
  // maximum number of faces per cell (to be found in proceeding loop)
  int nFacesPerCellMax = 0; 

  for (int i = 0; i < this->getNumberOfCells(); ++i)
  {
    // get cell i 
    dataSet->GetCell(i, genCell1);
    // get faces, find cells sharing it. if no cell shares it, 
    // use the locator of the surfWithPatches to find the patch number
    int numFaces = genCell1->GetNumberOfFaces();
    nFacesPerCellMax = (nFacesPerCellMax < numFaces ? numFaces : nFacesPerCellMax);
    for (int j = 0; j < numFaces; ++j)
    {
      vtkCell* face = genCell1->GetFace(j);
      bool shared = 0;
      int numVerts = face->GetNumberOfPoints();
      nVerticesPerFaceMax = (nVerticesPerFaceMax < numVerts ? numVerts : nVerticesPerFaceMax);
      facePtIds = face->GetPointIds(); 
      dataSet->GetCellNeighbors(i, facePtIds, sharedCellPtIds); 
      std::vector<int> facePntIds(numVerts);
      for (int k = 0; k < numVerts; ++k)
      {
        facePntIds[k] = face->GetPointId(k)+1;
      }
      if (sharedCellPtIds->GetNumberOfIds())
      {
        faceMap.insert(std::pair<std::vector<int>, std::pair<int,int>>
                        (facePntIds, std::make_pair(i+1, (int) sharedCellPtIds->GetId(0)+1))); 
      }
      else
      {
        double p1[3], p2[3], p3[3];
        face->GetPoints()->GetPoint(0,p1);
        face->GetPoints()->GetPoint(1,p2);
        face->GetPoints()->GetPoint(2,p3);
        double faceCenter[3];
        for (int k = 0; k < numVerts; ++k)
        {
          faceCenter[k] = (p1[k] + p2[k] + p3[k])/3.0;
        } 
        vtkIdType closestCellId;
        int subId;
        double minDist2;
        double closestPoint[3];
        // find closest point and closest cell to faceCenter
        surfCellLocator->FindClosestPoint(
          faceCenter, closestPoint, genCell2,closestCellId,subId,minDist2);
        double patchNo[1];
        surfWithPatches->getDataSet()->GetCellData()->GetArray("patchNo")
                                                  ->GetTuple(closestCellId, patchNo);
        faceMap.insert(std::pair<std::vector<int>, std::pair<int,int>>      
                  (facePntIds, std::make_pair(i+1, (int) -1*patchNo[0])));
      }
    }
  }
  
  std::map<int,int> patchMap;
  for (int i = 0; i < surfWithPatches->getNumberOfCells(); ++i)
  {
    double patchNo[1];
    surfWithPatches->getDataSet()->GetCellData()->GetArray("patchNo")
                                              ->GetTuple(i, patchNo);
    patchMap.insert(std::pair<int,int>(patchNo[0],i));
  }

  // write patch mapping file
  writePatchMap(mapFile, patchMap);
  // write cobalt file
  outputStream << 3 << "   " << 1 << "  " << patchMap.size() << std::endl;
  outputStream << this->getNumberOfPoints() << " " << faceMap.size()
               << " " << this->getNumberOfCells() << " "
               << nVerticesPerFaceMax << " " << nFacesPerCellMax << std::endl;
  for (int i = 0; i < this->getNumberOfPoints(); ++i)
  {
    std::vector<double> pnt(this->getPoint(i));
    outputStream << std::setw(21) << std::fixed << std::setprecision(15)
                 << pnt[0] << "   " << pnt[1] << "   " << pnt[2] << std::endl;
  }
  
  auto it = faceMap.begin();
  while (it != faceMap.end())
  {
    outputStream << it->first.size() << " ";
    auto faceIdIter = it->first.begin();
    while (faceIdIter != it->first.end())
    {
      outputStream << *faceIdIter << " ";
      ++faceIdIter;
    }
    outputStream << it->second.first << " " << it->second.second << std::endl;
    ++it;
  }

}

void meshBase::writeCobalt(meshBase* surfWithPatches, 
                           const std::string& mapFile, const std::string& ofname)
{
  std::ofstream outputStream(ofname);
  if(!outputStream.good()) 
  {
    std::cout << "Cannot open file " << ofname << std::endl;
    exit(1);
  }
  writeCobalt(surfWithPatches, mapFile,outputStream); 
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
              << " not found in set of point data arrays" << std::endl;
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
              << " not found in set of point data arrays" << std::endl;
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

bool sortIntVec_compare::operator() (std::vector<int> lhs, 
                                     std::vector<int> rhs) const
{
  std::sort(lhs.begin(), lhs.end());
  std::sort(rhs.begin(), rhs.end());
  return lhs < rhs;
}
