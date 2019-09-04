#include "ConversionDriver.H"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkModelMetadata.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

#include "AuxiliaryFunctions.H"
#include "cobalt.H"
#include "patran.H"
#include "pntMesh.H"
#include "vtkMesh.H"
#ifdef HAVE_EXODUSII
#  include "meshSrch.H"
#endif
#ifdef HAVE_CFMSH
#  include "foamMesh.H"
#  include "gmshMesh.H"
#endif

//----------------------- Conversion Driver
//-----------------------------------------//
ConversionDriver::ConversionDriver(const std::string &srcmsh,
                                   const std::string &trgmsh,
                                   const std::string &method,
                                   const std::string &ofname,
                                   const jsoncons::json &inputjson)
    : source(nullptr) {
  std::cout << "ConversionDriver created" << std::endl;
  nemAux::Timer T;
  T.start();
  // convert to pnt mesh
  if (method == "VTK->PNT") {
    source = meshBase::Create(srcmsh);
    // number of dimensions
    int dim = inputjson["Conversion Options"]["Dimension"].as<int>();
    // looping through blocks
    int nBlk = inputjson["Conversion Options"]["Block Data"].size();
    std::cout << "Number of Blocks : " << nBlk << std::endl;
    PNTMesh::BlockMap elmBlkMap;
    elmBlkMap.resize(nBlk);
    int iBlk = 0;
    for (auto it = inputjson["Conversion Options"]["Block Data"]
                       .array_range()
                       .begin();
         it !=
         inputjson["Conversion Options"]["Block Data"].array_range().end();
         ++it) {
      // element ids for block
      std::cout << (*it)["Name"].as<std::string>() << std::endl;
      std::vector<int> ids;
      if (it->contains("Element ID Range")) {
        int rs, re;
        rs = (*it)["Element ID Range"][0].as<int>();
        re = (*it)["Element ID Range"][1].as<int>();
        std::cout << "Range " << rs << " to " << re << std::endl;
        for (int eid = rs; eid <= re; eid++) ids.push_back(eid);
      } else {
        for (const auto &it2 : (*it)["Element ID"].array_range()) {
          ids.push_back(it2.as<int>());
          std::cout << it2.as<int>() << " ";
        }
        std::cout << std::endl;
      }
      // registering block information
      elmBlkMap[iBlk].ordIntrp = (*it)["Element Order"].as<int>();
      elmBlkMap[iBlk].ordEquat = (*it)["Equation Order"].as<int>();
      elmBlkMap[iBlk].eTpe =
          PNTMesh::elmTypeNum((*it)["Element Type"].as<std::string>());
      elmBlkMap[iBlk].regionName = (*it)["Name"].as<std::string>();
      std::string bcTagName = (*it)["BC Tag"].as<std::string>();
      elmBlkMap[iBlk].srfBCTag.push_back(PNTMesh::bcTagNum(bcTagName));
      elmBlkMap[iBlk].elmIds = ids;
      iBlk++;
    }
    auto *pm = new PNTMesh::pntMesh(source, dim, nBlk, elmBlkMap);
    pm->write(trgmsh);
    delete pm;
  } else if (method == "GMSH->EXO") {
#ifdef HAVE_EXODUSII  // NEMoSys is compiled with exodus
    // reading vitals
    std::cout << "Converting to EXODUS II..." << std::endl;
    jsoncons::json opts = inputjson["Conversion Options"];
    genExo(opts, ofname);
#else
    std::cerr << "Error: Compile NEMoSys with ENABLE_EXODUS to use this option."
              << std::endl;
    exit(-1);
#endif
  } else if (method == "FOAM->MSH") {
#ifdef HAVE_CFMSH
    meshBase *fm = new FOAM::foamMesh();
    fm->read("NULL");
    // TODO: Fix report and write methods for the foamMesh class
    // fm->setFileName(ofname);
    // fm->report();
    // fm->writeMSH();
    auto *gm = new gmshMesh(fm);
    gm->write(ofname);
    delete fm;
#else
    std::cerr << "Error: Compile NEMoSys with ENABLE_CFMSH to use this option."
              << std::endl;
    exit(-1);
#endif
  } else if (method == "FOAM->VTK") {
#ifdef HAVE_CFMSH
    // supports: vtu, vtk
    meshBase *fm = new FOAM::foamMesh();
    fm->read("NULL");
    // TODO: Fix report and write methods for the foamMesh class
    std::cout << "Variable values is = " << srcmsh << std::endl;
    vtkMesh *vm = new vtkMesh(fm->getDataSet(), ofname);
    vm->report();
    vm->write();
    delete vm;
    delete fm;
#else
    std::cerr << "Error: Compile NEMoSys with ENABLE_CFMSH to use this option."
              << std::endl;
    exit(-1);
#endif
  } else if (method == "GMSH->VTK") {
    if (srcmsh.find(".msh") != std::string::npos) {
      std::cout << "Detected file in GMSH format" << std::endl;
      std::cout << "Converting to VTK ...." << std::endl;
    } else {
      std::cerr << "Source mesh file is not in GMSH format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good()) {
      std::cerr << "Error opening file " << srcmsh << std::endl;
      exit(1);
    }

    std::string line;
    int numPoints, numCells;
    std::vector<std::vector<std::vector<double>>> pointData;
    std::vector<std::vector<std::vector<double>>> cellData;
    std::vector<std::string> pointDataNames;
    std::vector<std::string> cellDataNames;

    // declare points to be pushed into dataSet_tmp
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    // declare dataSet_tmp which will be associated to output vtkMesh
    vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    // map to hold true index of points (gmsh allows non-contiguous ordering)
    std::map<int, int> trueIndex;

    // declare array to store physical entity for each element
    vtkSmartPointer<vtkDataArray> physicalEntity =
        vtkSmartPointer<vtkIdTypeArray>::New();
    ;
    physicalEntity->SetNumberOfComponents(1);
    physicalEntity->SetName("patchNo");

    while (getline(meshStream, line)) {
      if (line.find("$Nodes") != std::string::npos) {
        getline(meshStream, line);
        std::stringstream ss(line);
        ss >> numPoints;
        int id;
        double x, y, z;
        // allocate memory for points
        points->SetNumberOfPoints(numPoints);
        for (int i = 0; i < numPoints; ++i) {
          getline(meshStream, line);
          std::stringstream ss(line);
          ss >> id >> x >> y >> z;
          double point[3];
          point[0] = x;
          point[1] = y;
          point[2] = z;
          // insert point i
          points->SetPoint(i, point);
          trueIndex.insert(std::pair<int, int>(id, i));
        }
        // inserting point array into dataSet_tmp
        dataSet_tmp->SetPoints(points);
      }

      if (line.find("$Elements") != std::string::npos) {
        getline(meshStream, line);
        std::stringstream ss(line);
        ss >> numCells;
        int id, type, numTags;
        // allocate space for cell connectivities
        dataSet_tmp->Allocate(numCells);
        for (int i = 0; i < numCells; ++i) {
          getline(meshStream, line);
          std::stringstream ss(line);
          ss >> id >> type >> numTags;
          vtkSmartPointer<vtkIdList> vtkcellIds =
              vtkSmartPointer<vtkIdList>::New();
          if (type == 2) {
            int tmp;
            double tmp2[1];
            for (int j = 0; j < numTags; ++j) {
              ss >> tmp2[0];
              if (j == 0) {
                physicalEntity->InsertNextTuple(tmp2);
              }
            }
            for (int j = 0; j < 3; ++j) {
              ss >> tmp;
              // insert connectivities for cell into cellIds container
              vtkcellIds->InsertNextId(trueIndex[tmp]);  //-1);
            }
            // insert connectivities for triangle elements into dataSet
            dataSet_tmp->InsertNextCell(VTK_TRIANGLE, vtkcellIds);
          } else if (type == 4) {
            int tmp;
            double tmp2[1];
            for (int j = 0; j < numTags; ++j) {
              ss >> tmp2[0];
              if (j == 0) {
                physicalEntity->InsertNextTuple(tmp2);
              }
            }
            for (int j = 0; j < 4; ++j) {
              ss >> tmp;
              // insert connectivities for cell into cellIds container
              vtkcellIds->InsertNextId(trueIndex[tmp]);  //-1);
            }
            // insert connectivities for tet elements into dataSet
            dataSet_tmp->InsertNextCell(VTK_TETRA, vtkcellIds);
          } else {
            std::cerr
                << "Only triangular and tetrahedral elements are supported!"
                << std::endl;
            exit(1);
          }
        }
      }

      if (line.find("$NodeData") != std::string::npos) {
        std::vector<std::vector<double>> currPointData;
        getline(meshStream, line);  // moving to num_string_tags
        std::stringstream ss(line);
        int num_string_tags;
        ss >> num_string_tags;
        std::string dataname;
        for (int i = 0; i < num_string_tags; ++i) {
          getline(meshStream, line);  // get string tag
          if (i == 0) {
            std::stringstream ss(line);
            ss >> dataname;
          }
        }
        dataname.erase(std::remove(dataname.begin(), dataname.end(), '"'),
                       dataname.end());
        pointDataNames.push_back(dataname);
        getline(meshStream, line);  // moving to num_real_tags
        std::stringstream ss1(line);
        int num_real_tags;
        ss1 >> num_real_tags;
        for (int i = 0; i < num_real_tags; ++i) getline(meshStream, line);

        getline(meshStream, line);  // moving to num_int_tags
        std::stringstream ss2(line);
        int num_int_tags;
        ss2 >> num_int_tags;
        int dt, dim = 0, numFields = 0, tmp;
        for (int i = 0; i < num_int_tags; ++i) {
          getline(meshStream, line);  // get int tag
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
        for (int i = 0; i < numFields; ++i) {
          std::vector<double> data(dim);
          int id;
          double val;
          getline(meshStream, line);
          std::stringstream ss(line);
          ss >> id;
          for (int j = 0; j < dim; ++j) {
            ss >> val;
            data[j] = val;
          }
          currPointData.push_back(data);
        }
        pointData.push_back(currPointData);
      }

      if (line.find("$ElementData") != std::string::npos) {
        std::vector<std::vector<double>> currCellData;
        getline(meshStream, line);  // moving to num_string_tags
        std::stringstream ss(line);
        int num_string_tags;
        ss >> num_string_tags;
        std::string dataname;
        for (int i = 0; i < num_string_tags; ++i) {
          getline(meshStream, line);  // get string tag
          if (i == 0) {
            std::stringstream ss(line);
            ss >> dataname;
          }
        }
        dataname.erase(std::remove(dataname.begin(), dataname.end(), '"'),
                       dataname.end());
        cellDataNames.push_back(dataname);
        getline(meshStream, line);  // moving to num_real_tags
        std::stringstream ss1(line);
        int num_real_tags;
        ss1 >> num_real_tags;
        for (int i = 0; i < num_real_tags; ++i) getline(meshStream, line);

        getline(meshStream, line);  // moving to num_int_tags
        std::stringstream ss2(line);
        int num_int_tags;
        ss2 >> num_int_tags;
        int dt, dim = 0, numFields = 0, tmp;
        for (int i = 0; i < num_int_tags; ++i) {
          getline(meshStream, line);  // get int tag
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
        for (int i = 0; i < numFields; ++i) {
          std::vector<double> data(dim);
          int id;
          double val;
          getline(meshStream, line);
          std::stringstream ss(line);
          ss >> id;
          for (int j = 0; j < dim; ++j) {
            ss >> val;
            data[j] = val;
          }
          currCellData.push_back(data);
        }
        cellData.push_back(currCellData);
      }
    }

    vtkMesh *outputMesh = new vtkMesh(dataSet_tmp, trgmsh);

    for (int i = 0; i < pointData.size(); ++i)
      outputMesh->setPointDataArray(&(pointDataNames[i])[0u], pointData[i]);
    for (int i = 0; i < cellData.size(); ++i)
      outputMesh->setCellDataArray(&(cellDataNames[i])[0u], cellData[i]);
    // add physical entity mesh to VTK
    outputMesh->getDataSet()->GetCellData()->AddArray(physicalEntity);
    outputMesh->write();
  } else if (method == "VTK->COBALT") {
    if (srcmsh.find(".vt") != std::string::npos) {
      std::cout << "Detected file in VTK format" << std::endl;
      std::cout << "Converting to COBALT ...." << std::endl;
    } else {
      std::cerr << "Source mesh file is not in VTK format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good()) {
      std::cerr << "Error opening file " << srcmsh << std::endl;
      exit(1);
    }

    // create meshBase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(srcmsh);

    // create Cobalt object from meshBase
    COBALT::cobalt *cm = new COBALT::cobalt(myMesh, srcmsh, trgmsh, ofname);
    // write to file
    cm->write();
  } else if (method == "VTK->PATRAN") {
    if (srcmsh.find(".vt") != std::string::npos) {
      std::cout << "Detected file in VTK format" << std::endl;
      std::cout << "Converting to PATRAN ...." << std::endl;
    } else {
      std::cerr << "Source mesh file is not in VTK format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good()) {
      std::cerr << "Error opening file " << srcmsh << std::endl;
      exit(1);
    }

    // looping through blocks
    int nBC = inputjson["Conversion Options"]["BC Info"].size();
    std::cout << "Number of Boundary Conditions read: " << nBC << std::endl;

    // PATRAN specific BC information
    // Each map below maps the indicated information from the "Patch Number"
    // specified in the file

    std::map<int, int> faceTypeMap;  // Rocfrac FSI Type;
    // 0 = no FSI, 1 = FSI w/ burn, 2 = FSI w/o burn, etc.
    // see Rocfrac manual for details
    std::map<int, int>
        nodeTypeMap;  // patch numbers as specified in RocfracControl.txt
    std::map<int, bool> nodeStructuralMap;  // boolean indicating structural BC
    std::map<int, bool> nodeMeshMotionMap;  // boolean indicating mesh motion BC
    std::map<int, bool> nodeThermalMap;  // boolean indicating heat transfer BC

    for (auto it =
             inputjson["Conversion Options"]["BC Info"].array_range().begin();
         it != inputjson["Conversion Options"]["BC Info"].array_range().end();
         ++it) {
      if ((*it)["BC Type"].as<std::string>() == "Face") {
        faceTypeMap[(*it)["Patch Number"].as<int>()] =
            (*it)["Rocfrac FSI Type"].as<int>();
      } else if ((*it)["BC Type"].as<std::string>() == "Node") {
        nodeTypeMap[(*it)["Patch Number"].as<int>()] =
            (*it)["RocfracControl Type"].as<int>();
        nodeStructuralMap[(*it)["Patch Number"].as<int>()] =
            (*it)["Structural"].as<bool>();
        nodeMeshMotionMap[(*it)["Patch Number"].as<int>()] =
            (*it)["Mesh Motion"].as<bool>();
        nodeThermalMap[(*it)["Patch Number"].as<int>()] =
            (*it)["Thermal"].as<bool>();
      }
    }

    // Accumulate node patch preference
    // (determines which patch a node belongs to if it borders two or more
    // patches)
    std::vector<int> nppVec;
    for (const auto &nppItr :
         inputjson["Conversion Options"]["Node Patch Preference"]
             .array_range()) {
      nppVec.push_back(nppItr.as<int>());
    }

    // create meshBase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(srcmsh);

    // create PATRAN object from meshBase
    auto *pat = new PATRAN::patran(myMesh, srcmsh, trgmsh, faceTypeMap,
                                   nodeTypeMap, nodeStructuralMap,
                                   nodeMeshMotionMap, nodeThermalMap, nppVec);
    // pat->write();
  } else if (method == "VTK->FOAM") {
  #ifdef HAVE_CFMSH
    if (srcmsh.find(".vt") != std::string::npos) {
      std::cout << "Detected file in VTK format" << std::endl;
      std::cout << "Converting to FOAM Mesh ...." << std::endl;
    } else {
      std::cerr << "Source mesh file is not in VTK format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good()) {
      std::cerr << "Error opening file " << srcmsh << std::endl;
      exit(1);
    }
    // create meshBase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(srcmsh);

    // create foamMesh object
    FOAM::foamMesh *fm = new FOAM::foamMesh(myMesh);

    // Write polyMesh
    // fm->report();
    fm->write(ofname);
  #else
    std::cerr << "Error: Compile NEMoSys with ENABLE_CFMSH to use this option."
              << std::endl;
    exit(-1);
  #endif
  }
  else if(method == "VTK-HEX->VTK-TET")
  {
    if (srcmsh.find(".vt") != std::string::npos)
    {
      std::cout << "Detected file in VTK format" << std::endl;
      std::cout << "Converting to tet Mesh ...." << std::endl;
    }
    else
    {
      std::cerr << "Source mesh file is not in VTK format" << std::endl;
    }

    std::ifstream meshStream(srcmsh);
    if (!meshStream.good())
    {
      std::cerr << "Error opening file " << srcmsh << std::endl;
      exit(1);
    }
    // create meshBase object
    std::shared_ptr<meshBase> myMesh = meshBase::CreateShared(srcmsh);

    // Converts hex mesh to tet mesh and writes in VTU file.
    myMesh->convertHexToTetVTK(myMesh->getDataSet());
    myMesh->report();
    myMesh->write(ofname);
  }
  else 
  {
    std::cerr << "Error: Conversion method " << method
              << " is not a valid option." << std::endl;
    exit(-1);
  }

  T.stop();
}

ConversionDriver::~ConversionDriver() {
  delete source;
  std::cout << "ConversionDriver destroyed" << std::endl;
}

ConversionDriver *ConversionDriver::readJSON(const jsoncons::json &inputjson) {
  std::cout << "Reading JSON object" << std::endl;

  std::string srcmsh;
  std::string trgmsh;
  std::string outmsh;
  std::string method;

  if (inputjson.contains("Mesh File Options")) {
    if (inputjson["Mesh File Options"].contains("Input Mesh Files")) {
      if (inputjson["Mesh File Options"]["Input Mesh Files"].contains(
              "Source Mesh"))
        srcmsh =
            inputjson["Mesh File Options"]["Input Mesh Files"]["Source Mesh"]
                .as<std::string>();

      if (inputjson["Mesh File Options"]["Input Mesh Files"].contains(
              "Target Mesh"))
        trgmsh =
            inputjson["Mesh File Options"]["Input Mesh Files"]["Target Mesh"]
                .as<std::string>();

      if (inputjson["Mesh File Options"].contains("Output Mesh File"))
        outmsh = inputjson["Mesh File Options"]["Output Mesh File"]
                     .as<std::string>();
    }
  }

  // minimal json input is "Conversion Options"->"Method"
  method = inputjson["Conversion Options"]["Method"].as<std::string>();

  // starting proper conversion driver, right now just one
  auto *convdrvobj =
      new ConversionDriver(srcmsh, trgmsh, method, outmsh, inputjson);
  return convdrvobj;
}

ConversionDriver *ConversionDriver::readJSON(const std::string &ifname) {
  std::cout << "Reading JSON file" << std::endl;

  std::ifstream inputStream(ifname);
  if (!inputStream.good() || nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  // checking if array
  if (inputjson.is_array()) {
    std::cerr
        << "Warning: Input is an array. Only first element will be processed"
        << std::endl;
    return ConversionDriver::readJSON(inputjson[0]);
  } else {
    return ConversionDriver::readJSON(inputjson);
  }
}

#ifdef HAVE_EXODUSII
void ConversionDriver::genExo(const jsoncons::json &opts,
                              const std::string &fname) {
  int nMsh = opts.get_with_default("Number of Mesh", 0);
  bool needsPP = opts.get_with_default("Post Processing", false);

  // sanity check
  if (nMsh == 0) {
    std::cerr << "Error: At least one mesh should be provided!" << std::endl;
    exit(-1);
  }

  // starting conversion operation
  auto *em = new NEM::MSH::EXOMesh::exoMesh(fname);

  // reading meshes
  int ieb = 0;          // Element Block counter.
  int ins = 0;          // Node Set counter.
  int iss = 0;          // Side Set counter.
  int ndeIdOffset = 0;  // for element blocks and node sets
  int elmIdOffset = 0;  // for side sets
  for (int iMsh = 0; iMsh < nMsh; iMsh++) {
    jsoncons::json mshOpts = opts["Mesh Data"];
    std::string mshFName = mshOpts[iMsh].get_with_default("File", "");
    std::string mshName = mshOpts[iMsh].get_with_default("Name", "default");
    bool usePhys = mshOpts[iMsh].get_with_default("Use Physical Groups", false);

    // reading input mesh
    std::string fileExt = nemAux::find_ext(mshFName);
    if (fileExt == ".g" || fileExt == ".e" || fileExt == ".exo" ||
        fileExt == ".exodus") {
      // if already exodus format, stitch and continue.
      NEM::MSH::EXOMesh::exoMesh em2;
      em2.read(mshFName);

      // Allow changing element block names using key-value.
      if (mshOpts[iMsh].contains("Element Block Names")) {
        std::map<std::string, std::string> elmBlkNames =
            mshOpts[iMsh]["Element Block Names"]
                .as<std::map<std::string, std::string>>();

        for (const auto &kv : elmBlkNames)
          em2.setBlockName(std::stoi(kv.first) - 1, kv.second);
      }

      // Allow defining a global node set containing all nodes.
      if (mshOpts[iMsh].contains("Add Global Node Set")) {
        NEM::MSH::EXOMesh::ndeSetType ns;
        ns.id = ++ins;
        ns.nNde = em2.getNumberOfNodes();
        ns.name = mshOpts[iMsh]["Add Global Node Set"].as<std::string>();
        ns.ndeIdOffset = ndeIdOffset;
        for (int iNde = 0; iNde < ns.nNde; ++iNde)
          ns.ndeIds.emplace_back(iNde + 1);
        em2.addNdeSet(ns);
      }

      em->stitch(em2);
      ieb += em2.getNumberOfElementBlocks();
      ins += em2.getNumberOfNodeSets();
      iss += em2.getNumberOfSideSets();
      ndeIdOffset += em2.getNumberOfNodes();
      elmIdOffset += em2.getNumberOfElements();
      continue;
    }

    meshBase *mb = meshBase::Create(mshFName);
    // mb->write("exo_inp_mesh_" + std::to_string(iMsh) + ".vtu");

    // adding information to exodusII object
    int ndeIdOffset_local = 0;
    int elmIdOffset_local = 0;

    // add nodes to database
    for (nemId_t iNde = 0; iNde < mb->getNumberOfPoints(); ++iNde)
      em->addNde(mb->getPoint(iNde));

    // node coordinate to one nodeSet
    NEM::MSH::EXOMesh::ndeSetType ns;
    ns.id = ++ins;
    ns.nNde = mb->getNumberOfPoints();
    ns.name = mshName;
    ns.ndeIdOffset = ndeIdOffset;
    for (int iNde = 0; iNde < ns.nNde; iNde++) {
      ns.ndeIds.emplace_back(iNde + 1);
    }
    em->addNdeSet(ns);
    ndeIdOffset_local += ns.nNde;

    // add element blocks

    // Element bucket where they will be sorted.
    // Outer layer: Specifies grouping, such as physical groups. The 0th index
    // is reserved for elements without a group.
    // Inner layer: Sorted by element type. The 0th index is reserved for
    // unsupported elements. The current support respects the ordering of the
    // NEM::MSH::EXOMesh::elementType enum:
    //     OTHER, TRIANGLE, QUAD, TETRA, HEX
    std::map<int, std::map<NEM::MSH::EXOMesh::elementType, std::vector<int>>>
        elmBucket;

    std::vector<double> grpIds(mb->getNumberOfCells(), 0.0);
    if (usePhys) mb->getCellDataArray("PhysGrpId", grpIds);

    for (int iElm = 0; iElm < mb->getNumberOfCells(); iElm++) {
      VTKCellType vtkType =
          static_cast<VTKCellType>(mb->getDataSet()->GetCellType(iElm));
      NEM::MSH::EXOMesh::elementType exoType =
          NEM::MSH::EXOMesh::v2eEMap(vtkType);
      elmBucket[grpIds[iElm]][exoType].emplace_back(iElm);
    }

    // sanity check
    int numUnsupported = 0;
    for (const auto &elmGroup : elmBucket)
      if (elmGroup.second.count(NEM::MSH::EXOMesh::elementType::WEDGE) != 0 ||
          elmGroup.second.count(NEM::MSH::EXOMesh::elementType::OTHER) != 0)
        numUnsupported +=
            elmGroup.second.at(NEM::MSH::EXOMesh::elementType::WEDGE).size() +
            elmGroup.second.at(NEM::MSH::EXOMesh::elementType::OTHER).size();
    if (numUnsupported > 0) {
      std::cerr << "WARNING: Detected " << numUnsupported
                << " unsupported elements.\n";
      throw;
    }

    // for each group and supported type, if existent, add an element block
    for (const auto &elmGroup : elmBucket) {
      for (const auto &elmIds : elmGroup.second) {
        if (elmIds.second.empty()) continue;  // skip if empty

        NEM::MSH::EXOMesh::elmBlkType eb;
        eb.id = ++ieb;
        eb.ndeIdOffset = ndeIdOffset;
        eb.nElm = elmIds.second.size();

        if (usePhys)
          eb.name = mshName + "_PhysGrp_" + std::to_string(ieb);
        else
          eb.name = mshName + "_" + std::to_string(ieb);

        switch (elmIds.first) {
          case NEM::MSH::EXOMesh::elementType::TRIANGLE:
            std::cout << "Number of triangular elements = " << eb.nElm << "\n";
            eb.ndePerElm = 3;
            eb.eTpe = NEM::MSH::EXOMesh::elementType::TRIANGLE;
            break;
          case NEM::MSH::EXOMesh::elementType::QUAD:
            std::cout << "Number of quadrilateral elements = " << eb.nElm
                      << "\n";
            eb.ndePerElm = 4;
            eb.eTpe = NEM::MSH::EXOMesh::elementType::QUAD;
            break;
          case NEM::MSH::EXOMesh::elementType::TETRA:
            std::cout << "Number of tetrahedral elements = " << eb.nElm << "\n";
            eb.ndePerElm = 4;
            eb.eTpe = NEM::MSH::EXOMesh::elementType::TETRA;
            break;
          case NEM::MSH::EXOMesh::elementType::HEX:
            std::cout << "Number of hexahedral elements = " << eb.nElm << "\n";
            eb.ndePerElm = 8;
            eb.eTpe = NEM::MSH::EXOMesh::elementType::HEX;
            break;
          case NEM::MSH::EXOMesh::elementType::WEDGE:
          case NEM::MSH::EXOMesh::elementType::OTHER:
          default:
            std::cerr << "WARNING: Processing unsupported element. Previous "
                         "sanity check failed!\n";
            throw;
        }

        eb.conn.reserve(eb.nElm * eb.ndePerElm);
        vtkIdList *nids = vtkIdList::New();
        for (const auto &iElm : elmIds.second) {
          mb->getDataSet()->GetCellPoints(iElm, nids);
          for (int in = 0; in < eb.ndePerElm; ++in) {
            // offset node ids by 1
            eb.conn.emplace_back(nids->GetId(in) + 1);
          }
        }

        std::cout << "Min node index = "
                  << *min_element(eb.conn.begin(), eb.conn.end()) << "\n"
                  << "Max node index = "
                  << *max_element(eb.conn.begin(), eb.conn.end()) << "\n";
        std::cout << "Starting node offset = " << ndeIdOffset << std::endl;

        em->addElmBlk(eb);
        elmIdOffset_local += eb.nElm;
      }
    }

    // add side sets
    if (mb->getMetadata()) {
      // get side set metadata
      vtkSmartPointer<vtkModelMetadata> metadata = mb->getMetadata();
      vtkSmartPointer<vtkStringArray> sdeSetNames = metadata->GetSideSetNames();
      int *sdeSetElmLst = metadata->GetSideSetElementList();
      int *sdeSetSdeLst = metadata->GetSideSetSideList();
      int *sdeSetSze = metadata->GetSideSetSize();

      for (int iSS = 0; iSS < metadata->GetNumberOfSideSets(); iSS++) {
        NEM::MSH::EXOMesh::sdeSetType ss;
        ss.id = ++iss;
        ss.name = sdeSetNames->GetValue(iSS);
        ss.nSde = sdeSetSze[iSS];
        ss.elmIds.assign(sdeSetElmLst, sdeSetElmLst + sdeSetSze[iSS]);
        ss.sdeIds.assign(sdeSetSdeLst, sdeSetSdeLst + sdeSetSze[iSS]);
        ss.elmIdOffset = elmIdOffset;
        em->addSdeSet(ss);

        // Advance pointer for reading next side set.
        sdeSetElmLst += sdeSetSze[iSS];
        sdeSetSdeLst += sdeSetSze[iSS];
      }
    }

    // offsetting starting ids for next file
    ndeIdOffset += ndeIdOffset_local;
    elmIdOffset += elmIdOffset_local;
    // clean up
    delete mb;
  }

  // writing the file
  em->write();
  em->report();

  // performing post-processing tasks
  if (needsPP) {
    int nTsk = opts.get_with_default("Number of Tasks", 0);
    jsoncons::json ppTsk = opts["Tasks"];
    for (int iTsk = 0; iTsk < nTsk; iTsk++) {
      std::string ppFName = ppTsk[iTsk].get_with_default("File", "");
      std::cout << "Reading Post Processing JSON file " << iTsk << std::endl;
      std::ifstream inputStream(ppFName);
      if (!inputStream.good() || nemAux::find_ext(ppFName) != ".json") {
        std::cerr << "Error opening file " << ppFName << std::endl;
        exit(1);
      }
      if (nemAux::find_ext(ppFName) != ".json") {
        std::cerr << "Input File must be in .json format" << std::endl;
        exit(1);
      }
      jsoncons::json ppJson;
      inputStream >> ppJson;
      procExo(ppJson, fname, em);
    }

    // writing augmented exo file
    em->write();
    em->report();
  }

  // clean up
  delete em;
}

void ConversionDriver::procExo(const jsoncons::json &ppJson,
                               const std::string &fname,
                               NEM::MSH::EXOMesh::exoMesh *em) {
  // converting to mesh base for geometric inquiry
  meshBase *mb = meshBase::Create(fname);

  // performing requested operation
  std::string opr = ppJson.get_with_default("Operation", "");
  if (opr == "Material Assignment") {
    meshSrch *ms = meshSrch::Create(mb);
    // gathering information about all zones
    // if densities are defined, materials with higher density will be
    // prioritized
    bool appDen = ppJson.get_with_default("Apply Density", false);

    std::map<std::pair<double, std::string>, std::set<int>> zoneGeom;

    // loop over all zones
    for (const auto &zone : ppJson["Zones"].array_range()) {
      // assuming first element is zone information keyed by zone name
      // that we do not care about yet
      std::string matName = zone[0].get_with_default("Material Name", "N/A");
      std::string shape = zone[0].get_with_default("Shape", "N/A");
      std::cout << "Processing zone: " << zone.object_range().begin()->key()
                << "  Material: " << matName << "  Shape: " << shape;

      double density = 1.0;  // Default density is 1.0
      if (appDen) {
        density = zone[0].get_with_default("Density", 1.0);
        std::cout << "  Density: " << density;
      }
      std::cout << std::endl;

      std::vector<nemId_t> lst;
      if (shape == "Box") {
        // Box shape. Requires 3-vector of Min and Max in x, y, and z, resp.
        std::vector<double> bb;
        bb.push_back(zone[0]["Params"]["Min"][0].as<double>());
        bb.push_back(zone[0]["Params"]["Max"][0].as<double>());
        bb.push_back(zone[0]["Params"]["Min"][1].as<double>());
        bb.push_back(zone[0]["Params"]["Max"][1].as<double>());
        bb.push_back(zone[0]["Params"]["Min"][2].as<double>());
        bb.push_back(zone[0]["Params"]["Max"][2].as<double>());

        ms->FindCellsWithinBounds(bb, lst, true);
      } else if (shape == "STL") {
        // STL shape. Only supports Tri surface.
        // Node Coordinates are given as 3-vector in an array.
        // Connectivities are given as 3-vectors in an array.
        // Tris are 0-indexed.
        std::vector<std::vector<double>> crds;
        std::vector<std::vector<vtkIdType>> conns;
        for (const auto &crd :
             zone[0]["Params"]["Node Coordinates"].array_range())
          crds.push_back(crd.as<std::vector<double>>());
        for (const auto &conn :
             zone[0]["Params"]["Connectivities"].array_range())
          conns.push_back(conn.as<std::vector<vtkIdType>>());

        ms->FindCellsInTriSrf(crds, conns, lst);
      } else if (shape == "Sphere") {
        // Sphere shape.
        // Center is a 3-vector in an array.
        // Radius is a double.
        std::vector<double> center =
            zone[0]["Params"]["Center"].as<std::vector<double>>();
        auto radius = zone[0]["Params"]["Radius"].as<double>();

        ms->FindCellsInSphere(center, radius, lst);
      } else {
        std::cerr << "WARNING: Skipping unknown zone shape: " << shape
                  << std::endl;
        continue;
      }

      if (zone[0].contains("Only From Block")) {
        std::string blkName = zone[0]["Only From Block"].as<std::string>();

        std::vector<std::string> elmBlkNames = em->getElmBlkNames();
        auto elmBlkName =
            std::find(elmBlkNames.begin(), elmBlkNames.end(), blkName);
        if (elmBlkName == elmBlkNames.end()) {
          std::cerr << "WARNING: Only From Block " << blkName
                    << " matches no available blocks. Continuing with no "
                       "restriction.\n";
        } else {
          std::vector<int> lst_int(lst.begin(), lst.end());
          bool allIn = false;
          lst_int = em->lstElmInBlk(
              std::distance(elmBlkNames.begin(), elmBlkName), lst_int, allIn);
          lst.assign(lst_int.begin(), lst_int.end());
        }
      }
      zoneGeom[{1.0 / density, matName}].insert(lst.begin(), lst.end());
    }

    if (appDen)
      std::cout << "Applying material zones based on density ordering"
                << std::endl;

    // adjusting exodus database accordingly
    for (const auto &zone : zoneGeom) {
      // zone is ((density, material name), ids)
      std::vector<int> elmLst;
      std::cout << "Manipulating ExodusDB for " << zone.first.second
                << std::endl;
      elmLst.insert(elmLst.end(), zone.second.begin(), zone.second.end());
      em->addElmBlkByElmIdLst(zone.first.second, elmLst);
    }
  } else if (opr == "Check Duplicate Elements") {
    std::cout << "Checking for existence of duplicate elements ... ";
    meshSrch *ms = meshSrch::Create(mb);
    bool ret = ms->chkDuplElm();
    if (ret) {
      std::cerr << " The exodus database contains duplicate elements."
                << std::endl;
      exit(-1);
    } else {
      std::cout << "False" << std::endl;
    }
  } else if (opr == "Remove Block") {
    std::string blkName = ppJson.get_with_default("Block Name", "");
    std::cout << "Removing Block " << blkName << std::endl;
    em->removeElmBlkByName(blkName);
  } else if (opr == "Snap Node Coords To Zero") {
    double tol = ppJson.get_with_default("Tolerance", 0.0);
    std::cout << "Snapping nodal coordinates to zero using tolerance: " << tol
              << std::endl;
    em->snapNdeCrdsZero(tol);
  } else if (opr == "Boundary Condition Assignment") {
    // For EP16 boundary conditions are simply translated to node sets. Node
    // sets may have shared nodes. In that case a node the order of nodeset
    // matter. A later node set supersedes an earlier one.
    meshSrch *ms = meshSrch::Create(mb);

    // gathering information about all boundary node sets
    jsoncons::json bcs = ppJson["Condition"];
    for (const auto &bc : bcs.array_range()) {
      std::set<nemId_t> pntIds;

      // identify node ids on each boundary
      std::string bcName = bc["Name"].as<std::string>();
      std::string bcTyp = bc["Boundary Type"].as<std::string>();
      if (bcTyp == "Faces") {
        std::vector<double> srfCrd;
        std::vector<nemId_t> srfConn;
        jsoncons::json nc = bc["Params"]["Node Coordinates"];
        for (const auto &crds : nc.array_range())
          for (const auto &cmp : crds.array_range())
            srfCrd.push_back(cmp.as<double>());
        jsoncons::json conn = bc["Params"]["Connectivities"];
        for (const auto &tri : conn.array_range())
          for (const auto &idx : tri.array_range())
            srfConn.push_back(idx.as<nemId_t>());
        ms->FindPntsOnTriSrf(srfCrd, srfConn, pntIds);
        std::cout << "Number of points residing on the boundary " << bcName
                  << " is " << pntIds.size() << std::endl;
      } else if (bcTyp == "Edges") {
        std::vector<double> edgeCrd;
        jsoncons::json ncs = bc["Params"]["Start"];
        for (const auto &crds : ncs.array_range())
          for (const auto &cmp : crds.array_range())
            edgeCrd.push_back(cmp.as<double>());
        jsoncons::json nce = bc["Params"]["End"];
        for (const auto &crds : nce.array_range())
          for (const auto &cmp : crds.array_range())
            edgeCrd.push_back(cmp.as<double>());
        ms->FindPntsOnEdge(edgeCrd, pntIds);
        std::cout << "Number of points residing on the boundary " << bcName
                  << " is " << pntIds.size() << std::endl;
      } else {
        std::cerr << "Warning: unsupported boundary type " << bcTyp
                  << std::endl;
      }

      // register node set in Exodus II database
      if (!pntIds.empty()) {
        std::vector<int> nv;
        std::copy(pntIds.begin(), pntIds.end(), std::back_inserter(nv));
        em->addNdeSetByNdeIdLst(bcName, nv);
      }
    }
  } else if (opr == "Merge Nodes") {
    em->mergeNodes();
  } else {
    std::cerr << "Unknown operation requested: " << opr << std::endl;
  }
}

#endif  // HAVE_EXODUSII
