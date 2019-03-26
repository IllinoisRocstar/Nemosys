#include <algorithm>
#include <map>
#include <unordered_map>

// Nemosys
#include <ConversionDriver.H>
#include "meshSrch.H"
#include <AuxiliaryFunctions.H>
#include <vtkMesh.H>
#include <vtkCellData.h>
#include <vtkIdList.h>


// vtk
#include <vtkIdList.h>
#include <vtkCellData.h>

//----------------------- Conversion Driver -----------------------------------------//
ConversionDriver::ConversionDriver(std::string srcmsh, std::string trgmsh,
                               std::string method, std::string ofname,
                               json inputjson)
    : source(NULL)
{
  std::cout << "ConversionDriver created" << std::endl;
  Timer T;
  T.start();
  // convert to pnt mesh
  if (!method.compare("VTK->PNT"))
  {
    source = meshBase::Create(srcmsh);
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
    delete pm;
    pm = NULL;
  }
  else if (!method.compare("EXODUSII"))
  {
#ifdef HAVE_EXODUSII // NEMOSYS is compiled with exodus
    // reading vitals
    std::cout << "Converting to EXODUS II...\n";
    json opts = inputjson["Conversion Options"];
    int nMsh = opts.get_with_default("Number of Mesh", 0);
    bool needsPP = opts.get_with_default("Post Processing", false);
    
    // sanity check
    if (nMsh == 0)
    {
        std::cerr << "Error: At least one mesh should be provided!\n";
        exit(-1);
    }

    // starting conversion operation
    EXOMesh::exoMesh* em = new EXOMesh::exoMesh(ofname);
    
    // reading meshes
    int ndeIdOffset = 0;
    int ins = 0;
    int ieb = 0;
    for (int iMsh=0; iMsh < nMsh; iMsh++)
    {
        json mshOpts = opts["Mesh Data"];
        std::string mshFName = mshOpts[iMsh].get_with_default("File", "");
        std::string mshName = mshOpts[iMsh].get_with_default("Name", "default");
        bool usePhys = mshOpts[iMsh].get_with_default("Use Physical Groups", false);

        // reading input mesh
        meshBase* mb = meshBase::Create(mshFName);
        //mb->write("exo_inp_mesh_"+std::to_string(iMsh)+".vtu");

        // adding information to exodusII object
        if (!usePhys)
        {
            // node coordinate to one nodeSet
            EXOMesh::ndeSetType ns;
            ns.id = ++ins;
            ns.name = mshName;
            ns.usrNdeIds = false; // automatically node ids
            ns.nNde = mb->getNumberOfPoints();
            for (int iNde=0; iNde<ns.nNde; iNde++)
                ns.crds.push_back(mb->getPoint(iNde));
            em->addNdeSet(ns);

            // all elements into one element block
            // assuming input mesh has only one element type
            EXOMesh::elmBlockType eb;
            eb.id = ++ieb;
            eb.name = mshName;
            eb.ndeIdOffset = ndeIdOffset;
            int iTet = 0; 
            int iHex = 0; 
            for (int iElm=0; iElm<mb->getNumberOfCells(); iElm++)
            {
                EXOMesh::elementType eTpe = EXOMesh::v2eEMap((VTKCellType) (mb->getDataSet()->GetCellType(iElm)));
                if (eTpe == EXOMesh::elementType::TETRA)
                {
                    iTet++;
                }
                else if (eTpe == EXOMesh::elementType::HEX)
                {
                    iHex++;
                }
                else
                    continue;
                vtkIdList* nids = vtkIdList::New();
                mb->getDataSet()->GetCellPoints(iElm, nids);
                for (int in=0; in< nids->GetNumberOfIds(); in++)
                    eb.conn.push_back(nids->GetId(in)+1);  // offset node ids by 1
            }
            // sanity check
            if (iTet*iHex != 0)
            {
                std::cerr << "Mixed Hex/Tet meshes are not supported!\n";
                throw;
            }
            if (iTet != 0)
                std::cout << "Number of tetrahedral elements = " << iTet << std::endl;
            else
                std::cout << "Number of hexahedral elements = " << iHex << std::endl;
            std::cout << "Min nde indx = " << *min_element(eb.conn.begin(), eb.conn.end()) << "\n"
                      << "Max nde indx = " << *max_element(eb.conn.begin(), eb.conn.end()) << "\n";
            std::cout << "Starting node offset = " << ndeIdOffset << std::endl;
            if (iTet != 0)
            {
                eb.nElm = iTet;
                eb.eTpe = EXOMesh::elementType::TETRA;
                eb.ndePerElm = 4;
            } 
            else 
            {
                eb.nElm = iHex;
                eb.eTpe = EXOMesh::elementType::HEX;
                eb.ndePerElm = 8;
            }
            em->addElmBlk(eb);
            // offseting starting node id for next file
            ndeIdOffset += ns.nNde;
        }
        else
        {
            // get number of physical groups
            // loop through cell data and identify physical groups
            // we also filter only TETRA/HEX cells
            int nc = mb->getNumberOfCells();
            // check physical group (PhysGrpId) exists, obtain id, and data array
            int physGrpArrIdx = mb->getCellDataIdx("PhysGrpId");
            if (physGrpArrIdx < 0)
            {
                std::cerr << "Error : Input dataset does not have PhyGrpId cell data! Aborting.\n";
                exit(-1);
            }
            vtkCellData *cd = mb->getDataSet()->GetCellData();
            vtkDataArray* physGrpIds = cd->GetArray(physGrpArrIdx);
            
            // loop through elements and obtain physical group ids
            std::vector<int> elmIds;
            std::vector<int> elmPhysGrp;
            for (int ic=0; ic<nc; ic++)
            {
                EXOMesh::elementType eTpe = EXOMesh::v2eEMap((VTKCellType) (mb->getDataSet()->GetCellType(ic)));
                if (eTpe != EXOMesh::elementType::TETRA && eTpe != EXOMesh::elementType::HEX)
                    continue;
                elmIds.push_back(ic);
                double* tmp = physGrpIds->GetTuple(ic);
                elmPhysGrp.push_back(int(*tmp));
            }

            // number of unique physical groups
            std::set<int> unqPhysGrpIds;
            int nPhysGrp;
            for (auto it = elmPhysGrp.begin(); it != elmPhysGrp.end(); it++)
                unqPhysGrpIds.insert(*it);
            nPhysGrp = unqPhysGrpIds.size();
            std::cout << "Number of physical groups : " << nPhysGrp << std::endl;

            // one node set for all groups
            // node coordinate to nodeSet
            EXOMesh::ndeSetType ns;
            ns.id = ++ins;
            ns.name = mshName;
            ns.usrNdeIds = false; // automatically node ids
            ns.nNde = mb->getNumberOfPoints();
            for (int iNde=0; iNde<ns.nNde; iNde++)
                ns.crds.push_back(mb->getPoint(iNde));
            em->addNdeSet(ns);
            
            // for each physical group one element block
            for (auto it1=unqPhysGrpIds.begin(); it1!=unqPhysGrpIds.end(); it1++)
            {
                // element to elementBlock
                EXOMesh::elmBlockType eb;
                eb.id = ++ieb;
                eb.name = mshName+"_PhysGrp_"+std::to_string(ieb);
                eb.ndeIdOffset = ndeIdOffset;
                int iTet = 0; 
                int iHex = 0;
                for (auto it2=elmIds.begin(); it2!=elmIds.end(); it2++)
                {
                    double* tmp2 = physGrpIds->GetTuple(*it2);
                    if (int(*tmp2) != *it1)
                        continue;
                    vtkIdList* nids = vtkIdList::New();
                    mb->getDataSet()->GetCellPoints(*it2, nids);
                    if (nids->GetNumberOfIds() == 4)
                        iTet++;
                    else if (nids->GetNumberOfIds() == 8)
                        iHex++;
                    for (int in=0; in< nids->GetNumberOfIds(); in++)
                        eb.conn.push_back(nids->GetId(in)+1);  // offset node ids by 1
                }
                // sanity check
                if (iTet*iHex != 0)
                {
                    std::cerr << "Mixed Hex/Tet meshes are not supported!\n";
                    throw;
                }
                if (iTet != 0)
                    std::cout << "Number of group tetrahedral elements = " << iTet << std::endl;
                else if (iHex != 0)
                    std::cout << "Number of group hexahedral elements = " << iHex << std::endl;
                std::cout << "Min nde indx = " << *min_element(eb.conn.begin(), eb.conn.end()) << "\n"
                          << "Max nde indx = " << *max_element(eb.conn.begin(), eb.conn.end()) << "\n";
                std::cout << "Starting node offset = " << ndeIdOffset << std::endl;
                if (iTet != 0)
                {
                    eb.nElm = iTet;
                    eb.eTpe = EXOMesh::elementType::TETRA;
                    eb.ndePerElm = 4;
                    em->addElmBlk(eb);
                }
                else
                {
                    eb.nElm = iHex;
                    eb.eTpe = EXOMesh::elementType::HEX;
                    eb.ndePerElm = 8;
                    em->addElmBlk(eb);
                }
            }
            // offset starting node id for next file
            ndeIdOffset += ns.nNde;
        }
        // clean up
        delete mb;
        mb = NULL;
    }

    // writing the file
    em->write();
    em->report();

    // performing post-processing tasks
    if (needsPP)
    {   
      int nTsk = opts.get_with_default("Number of Tasks", 0);
      json ppTsk = opts["Tasks"];
      for (int iTsk=0; iTsk<nTsk; iTsk++)
      {
        std::string ppFName= ppTsk[iTsk].get_with_default("File", "");
        std::cout << "Reading Post Processing JSON file "<< iTsk << std::endl;
        std::ifstream inputStream(ppFName);
        if (!inputStream.good() || find_ext(ppFName) != ".json")
        {
          std::cout << "Error opening file " << ppFName << std::endl;
          exit(1);
        }
        if (find_ext(ppFName) != ".json")
        {
          std::cout << "Input File must be in .json format" << std::endl;
          exit(1);
        }
        json ppJson;
        inputStream >> ppJson;
        procExo(ppJson, ofname, em);
      }

      // writing augmented exo file
      em->write();
      em->report();
    }
    
    // clean up
    delete em;
    em = NULL;
#else
    std::cerr << "Error: Compile NEMOSYS with ENABLE_EXODUS to use this option.\n";
    exit(-1);
#endif
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
 
  srcmsh = "";
  trgmsh = "";
  outmsh = "";

  if (inputjson.has_key("Mesh File Options"))
  {
      if (inputjson["Mesh File Options"].has_key("Input Mesh Files"))
      {
          if (inputjson["Mesh File Options"]["Input Mesh Files"].has_key("Source Mesh"))
            srcmsh = inputjson["Mesh File Options"]
                                ["Input Mesh Files"]
                                ["Source Mesh"].as<std::string>();

          if (inputjson["Mesh File Options"]["Input Mesh Files"].has_key("Target Mesh"))
            trgmsh = inputjson["Mesh File Options"]
                                ["Input Mesh Files"]
                                ["Target Mesh"].as<std::string>();
          
          if (inputjson["Mesh File Options"].has_key("Output Mesh File"))
            outmsh = inputjson["Mesh File Options"]
                                ["Output Mesh File"].as<std::string>();
      }
  }

  // minimal json input is "Conversion Options"->"Method"
  method = inputjson["Conversion Options"]
                    ["Method"].as<std::string>(); 

  // starting proper conversion driver, right now just one
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

#ifdef HAVE_EXODUSII
void ConversionDriver::procExo(json ppJson, std::string fname, EXOMesh::exoMesh* em)
{
  // converting to mesh base for geometric inquiry
  meshBase* mb = meshBase::Create(fname); 

  // performing requested operation
  std::string opr = ppJson.get_with_default("Operation", "");
  if (!opr.compare("Material Assignment"))
  {
      meshSrch* ms = meshSrch::Create(mb);
      // gathering information about all zones
      // if densities are defined materials with higher
      // density will be prioratized 
      bool appDen = ppJson.get_with_default("Apply Density", false);
      std::unordered_map<std::string, std::set<int> > zoneGeom;
      json zones = ppJson["Zones"];
      int nZn = zones.size();
      // ordering zones based on density 
      // default is 1.0
      std::map<double,std::vector<int> > znOrd;
      if (appDen)
      {
          std::cout << "Applying material zones based on density ordering\n";
          // order based of density
          for (int iZn=0; iZn<nZn; iZn++)
          {
              json znInfo = zones[iZn][0];
              double density = znInfo.get_with_default("Density",1.0);
              znOrd[density].push_back(iZn);
          }
      } 
      else 
      {
          // order based on apperance in the file
          for (int iZn=0; iZn<nZn; iZn++)
              znOrd[1.0].push_back(iZn);
      }

      for (auto im=znOrd.begin(); im!=znOrd.end(); im++)
      {
          for (auto iz=(im->second).begin(); iz!=(im->second).end(); iz++)
          {
              int iZn = *iz;
              // assuming first element is zone infomration
              // keyed by zone name that we do not care about
              // yet
              json znInfo = zones[iZn][0];
              std::string matName = znInfo.get_with_default("Material Name","N/A");
              std::string shape = znInfo.get_with_default("Shape","N/A");
              std::cout <<"Processing zone "<<iZn
                  <<" Material "<<matName
                  <<" Shape "<<shape<<std::endl;

              if (!shape.compare("Box"))
              {
                  std::vector<double> bb;
                  bb.push_back( znInfo["Params"]["Min"][0].as<double>() ); 
                  bb.push_back( znInfo["Params"]["Max"][0].as<double>() ); 
                  bb.push_back( znInfo["Params"]["Min"][1].as<double>() ); 
                  bb.push_back( znInfo["Params"]["Max"][1].as<double>() ); 
                  bb.push_back( znInfo["Params"]["Min"][2].as<double>() ); 
                  bb.push_back( znInfo["Params"]["Max"][2].as<double>() ); 
                  
                  std::vector<int> lst;
                  ms->FindCellsWithinBounds(bb, lst, true);
                  zoneGeom[matName].insert(lst.begin(), lst.end());
                  
              }
              else
                  std::cout << "WARNNING: Skipiing unkown zone shape: " 
                      << shape << std::endl;     
          }
      }

      // adjusting exodus database accordingly
      // unordered_maps act LIFO
      for (auto it1=zoneGeom.begin(); it1!=zoneGeom.end(); it1++)
      {
          std::vector<int> elmLst;
          std::cout << "Manipulating ExodusDB for " << (it1->first) << std::endl;
          elmLst.insert(elmLst.end(), (it1->second).begin(), (it1->second).end());
          em->addElmBlkByElmIdLst(it1->first, elmLst);
      }
  }
  else if (!opr.compare("Check Duplicate Elements"))
  {
      std::cout << "Checking for existance of duplicate elements ... ";
      meshSrch* ms = meshSrch::Create(mb);
      bool ret = ms->chkDuplElm(); 
      if (ret)
      {
          std::cerr << " The exodus database contains duplicate elements.\n";
          exit(-1);
      } 
      else 
        std::cout << "False\n";
  }
  else if (!opr.compare("Remove Block"))
  {
      std::string blkName = ppJson.get_with_default("Block Name", "");
      std::cout << "Removing Block " << blkName << std::endl;
      em->removeElmBlkByName(blkName);
  }
  else if (!opr.compare("Snap Node Coords To Zero"))
  {
      double tol = ppJson.get_with_default("Tolerance", 0.0);
      std::cout << "Snapping nodal coordinates to zero using tolerance " << tol << std::endl;
      em->snapNdeCrdsZero(tol);
  }
  else if (!opr.compare("Boundary Condition Assignment"))
  {
      // For EP16 boundary conditions are simply translated
      // to node sets. Node sets may have shared nodes. In
      // that case a node the order of nodeset matter. A later
      // node set superceeds an earlier one.
      meshSrch* ms = meshSrch::Create(mb);      
      // gathering information about all boundary node sets
      json bcs = ppJson["Condition"];
      int nBC = bcs.size();
      for (auto bc : bcs.array_range())
      {
        std::set<int> pntIds;              
        // identify node ids on each boundary
        std::string bcName = bc["Name"].as<std::string>();
        std::string bcTyp = bc["Boundary Type"].as<std::string>();
        if (!bcTyp.compare("Faces"))
        {
           std::vector<double> srfCrd;
           std::vector<int> srfConn;
           json nc = bc["Params"]["Node Coordinates"];
           for (auto crds : nc.array_range())
               for (auto cmp : crds.array_range())
                   srfCrd.push_back(cmp.as<double>());
           json conn = bc["Params"]["Connectivities"];
           for (auto tri : conn.array_range())
               for (auto idx : tri.array_range())
                   srfConn.push_back(idx.as<int>());
           ms->FindPntsOnTriSrf(srfCrd, srfConn, pntIds);
           std::cout << "Number of points resinding on the boundary " 
               << bcName << " is " << pntIds.size() << "\n";
        }
        else if (!bcTyp.compare("Edges"))
        {
           std::vector<double> edgeCrd;
           json ncs = bc["Params"]["Start"];
           for (auto crds : ncs.array_range())
               for (auto cmp : crds.array_range())
                   edgeCrd.push_back(cmp.as<double>());
           json nce = bc["Params"]["End"];
           for (auto crds : nce.array_range())
               for (auto cmp : crds.array_range())
                   edgeCrd.push_back(cmp.as<double>());
           ms->FindPntsOnEdge(edgeCrd, pntIds);
           std::cout << "Number of points resinding on the boundary " 
               << bcName << " is " << pntIds.size() << "\n";
        }
        else
            std::cerr << "Warning: unsupported bounday type " << bcTyp << std::endl;
        // register node set in Exodus II database
        if (!pntIds.empty())
        {
            std::vector<int> nv;
            std::copy(pntIds.begin(), pntIds.end(), std::back_inserter(nv));
            em->addNdeSetByNdeIdLst(bcName, nv);
        }
      }
  }
  else if (!opr.compare("Merge Nodes"))
  {
      em->mergeNodes();
  }
  else
  {
      std::cout << "Unknown operation requested : " << opr << std::endl;
  }
}
#endif

