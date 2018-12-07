#include <ConversionDriver.H>
#include <AuxiliaryFunctions.H>
#include <vtkMesh.H>
#include <vtkCellData.h>


//----------------------- Conversion Driver -----------------------------------------//
ConversionDriver::ConversionDriver(std::string srcmsh, std::string trgmsh,
                               std::string method, std::string ofname,
                               json inputjson)
{
  source = meshBase::Create(srcmsh);
  std::cout << "ConversionDriver created" << std::endl;
  Timer T;
  T.start();
  // convert to pnt mesh
  if (!method.compare("VTK->PNT"))
  {
    // number of dimensions
    int dim = inputjson["Conversion Options"]
                       ["Dimension"].as<int>();
    // looping through blocks
    int nBlk = inputjson["Conversion Options"]["Block Data"].size();
    std::cout << "Number of Blocks : " << nBlk << std::endl;
    PNTMesh::BlockMap elmBlkMap;
    elmBlkMap.resize(nBlk);
    int iBlk = 0;
    for (auto it = inputjson["Conversion Options"]
                            ["Block Data"].array_range().begin();
              it!= inputjson["Conversion Options"]
                            ["Block Data"].array_range().end();
              ++it)
    {
      // element ids for block
      std::cout << (*it)["Name"].as<std::string>() << std::endl;
      std::vector<int> ids;
      if (it->has_key("Element ID Range"))
      {
        int rs,re;
        rs = (*it)["Element ID Range"][0].as<int>();
        re = (*it)["Element ID Range"][1].as<int>();
        std::cout << "Range " << rs << " to " << re << std::endl;
        for (int eid= rs;
                 eid <= re;
                 eid++)
          ids.push_back(eid);
      }
      else
      {
        for (auto it2 = (*it)["Element ID"].array_range().begin();
                  it2 != (*it)["Element ID"].array_range().end();
                  it2++)
        {
          ids.push_back( (*it2).as<int>() );
          std::cout << (*it2).as<int>() << " ";
        }
        std::cout << std::endl;
      }
      // registering block information
      elmBlkMap[iBlk].ordIntrp = (*it)["Element Order"].as<int>();
      elmBlkMap[iBlk].ordEquat = (*it)["Equation Order"].as<int>();
      elmBlkMap[iBlk].eTpe = PNTMesh::elmTypeNum( (*it)["Element Type"].as<std::string>() );
      elmBlkMap[iBlk].regionName = (*it)["Name"].as<std::string>();
      std::string bcTagName = (*it)["BC Tag"].as<std::string>();
      elmBlkMap[iBlk].srfBCTag.push_back( PNTMesh::bcTagNum( bcTagName ) );
      elmBlkMap[iBlk].elmIds = ids;
      iBlk++;
    }
    PNTMesh::pntMesh* pm = new PNTMesh::pntMesh(source, dim, nBlk, elmBlkMap);
    pm->write(trgmsh);
  }
  else if (!method.compare("GMSH->VTK"))
  {
    if (srcmsh.find(".msh") != -1)
    {
      std::cout << "Detected file in GMSH format" << std::endl;
      std::cout << "Converting to VTK ...." << std::endl;
    }
    else
    {
      std::cout << "Source mesh file is not in GMSH format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good())
    {
      std::cout << "Error opening file " << srcmsh << std::endl;
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
  
    // declare array to store physical entity for each element
    vtkSmartPointer<vtkDataArray> physicalEntity = vtkSmartPointer<vtkIdTypeArray>::New();;
    physicalEntity->SetNumberOfComponents(1);
    physicalEntity->SetName("patchNo");
  
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
            double tmp2[1];
            for (int j = 0; j < numTags; ++j)
            {
              ss >> tmp2[0];
              if (j == 0)
              {
                physicalEntity->InsertNextTuple(tmp2);
              }
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
            double tmp2[1];
            for (int j = 0; j < numTags; ++j)
            {
              ss >> tmp2[0];
              if (j == 0)
              {
                physicalEntity->InsertNextTuple(tmp2);
              }
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
  
    vtkMesh* outputMesh = new vtkMesh(dataSet_tmp, trgmsh);

    for (int i = 0; i < pointData.size(); ++i)
      outputMesh->setPointDataArray(&(pointDataNames[i])[0u], pointData[i]);
    for (int i = 0; i < cellData.size(); ++i)
      outputMesh->setCellDataArray(&(cellDataNames[i])[0u], cellData[i]);
    // add physical entity mesh to VTK
    outputMesh->getDataSet()->GetCellData()->AddArray(physicalEntity);
    outputMesh->write();
  }
  else if (!method.compare("VTK->COBALT"))
  {
    if (srcmsh.find(".vt") != -1)
    {
      std::cout << "Detected file in VTK format" << std::endl;
      std::cout << "Converting to COBALT ...." << std::endl;
    }
    else
    {
      std::cout << "Source mesh file is not in VTK format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good())
    {
      std::cout << "Error opening file " << srcmsh << std::endl;
      exit(1);
    }  
    
    // create meshbase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(srcmsh);

    // create Cobalt object from meshBase
    COBALT::cobalt* cm = new COBALT::cobalt(myMesh, srcmsh, trgmsh, ofname);
    // write to file
    cm->write();
  }
  else if (!method.compare("VTK->PATRAN"))
  {
    if (srcmsh.find(".vt") != -1)
    {
      std::cout << "Detected file in VTK format" << std::endl;
      std::cout << "Converting to PATRAN ...." << std::endl;
    }
    else
    {
      std::cout << "Source mesh file is not in VTK format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good())
    {
      std::cout << "Error opening file " << srcmsh << std::endl;
      exit(1);
    }  
    
    // looping through blocks
    int nBC = inputjson["Conversion Options"]["BC Info"].size();
    std::cout << "Number of Boundary Conditions read: " << nBC << std::endl;

    // Patran specific BC information
    // Each map below maps the indicated information from the "Patch Number" specified in the file

    std::map<int, int> faceTypeMap;         // Rocfrac FSI Type;
                                            // 0 = no FSI, 1 = FSI w/ burn, 2 = FSI w/o burn, etc.
                                            // see Rocfrac manual for details
    std::map<int, int> nodeTypeMap;         // patch numbers as specified in RocfracControl.txt
    std::map<int, bool> nodeStructuralMap;  // boolean indicating structural BC
    std::map<int, bool> nodeMeshMotionMap;  // boolean indicating mesh motion BC
    std::map<int, bool> nodeThermalMap;     // boolean indicating heat transfer BC


    for (auto it = inputjson["Conversion Options"]
                            ["BC Info"].array_range().begin();
              it!= inputjson["Conversion Options"]
                            ["BC Info"].array_range().end();
              ++it)
    {
      if ((*it)["BC Type"].as<std::string>() == "Face")
      {
        faceTypeMap[(*it)["Patch Number"].as<int>()] = (*it)["Rocfrac FSI Type"].as<int>();
      }
      else if ((*it)["BC Type"].as<std::string>() == "Node")
      {
        nodeTypeMap[(*it)["Patch Number"].as<int>()] = (*it)["RocfracControl Type"].as<int>();
        nodeStructuralMap[(*it)["Patch Number"].as<int>()] = (*it)["Structural"].as<bool>();
        nodeMeshMotionMap[(*it)["Patch Number"].as<int>()] = (*it)["Mesh Motion"].as<bool>();
        nodeThermalMap[(*it)["Patch Number"].as<int>()] = (*it)["Thermal"].as<bool>();
      }
    }
    
    // Accumulate node patch preference (determines which patch a node belongs to if it borders two or more patches)
    std::vector<int> nppVec;
    for (auto nppItr = inputjson["Conversion Options"]["Node Patch Preference"].array_range().begin();
      nppItr != inputjson["Conversion Options"]["Node Patch Preference"].array_range().end(); ++nppItr)
    {
      nppVec.push_back((*nppItr).as<int>());
    }

    // create meshbase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(srcmsh);

    // create Patran object from meshBase
    PATRAN::patran* pat = new PATRAN::patran(myMesh, srcmsh, trgmsh,
                                            faceTypeMap, nodeTypeMap,
                                            nodeStructuralMap, nodeMeshMotionMap,
                                            nodeThermalMap, nppVec);
    //pat->write();
  }

  T.stop();
}


ConversionDriver::~ConversionDriver()
{
  if (source)
  {
    delete source;
    source = 0;
  }
  std::cout << "ConversionDriver destroyed" << std::endl;
}

ConversionDriver* ConversionDriver::readJSON(json inputjson)
{

  std::cout << "Reading JSON object\n";
  std::string srcmsh; 
  std::string trgmsh; 
  std::string outmsh; 
  std::string method;
  
  srcmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Source Mesh"].as<std::string>();
  trgmsh = inputjson["Mesh File Options"]
                    ["Input Mesh Files"]
                    ["Target Mesh"].as<std::string>();
  outmsh = inputjson["Mesh File Options"]
                    ["Output Mesh File"].as<std::string>();
  method = inputjson["Conversion Options"]
                    ["Method"].as<std::string>(); 

  ConversionDriver* convdrvobj;
  convdrvobj = new ConversionDriver(srcmsh, trgmsh, method, outmsh, inputjson);
  return convdrvobj;

}

ConversionDriver* ConversionDriver::readJSON(std::string ifname)
{
  std::cout << "Reading JSON file\n";
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || find_ext(ifname) != ".json")
  {
    std::cout << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (find_ext(ifname) != ".json")
  {
    std::cout << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  json inputjson;
  inputStream >> inputjson;
  
  // checking if array
  if (inputjson.is_array())
  {
    std::cout << "Warning: Input is an array. Only first element will be processed\n";
    return ConversionDriver::readJSON(inputjson[0]);
  } 
  else
  {
    return ConversionDriver::readJSON(inputjson);
  }
    
}
