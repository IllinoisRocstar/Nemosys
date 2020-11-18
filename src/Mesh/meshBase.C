#include "meshBase.H"

#include <exodusII.h>
#include <vtkAppendFilter.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
#include <vtkDataArray.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkExtractSelection.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkUnstructuredGrid.h>

//#include "cobalt.H"
//#include "patran.H"
#include "AuxiliaryFunctions.H"
#include "Cubature.H"
#include "MeshQuality.H"
#include "Refine.H"
#include "SizeFieldGen.H"
#include "exoMesh.H"
#include "meshGen.H"
#include "meshPartitioner.H"
#include "pntMesh.H"
#include "vtkMesh.H"

// netgen
#ifdef HAVE_NGEN
namespace nglib {
#  include <nglib.h>
}
#endif

// TODO: Stop using setPoint/CellDataArray in export methods
//        - instead, use the faster vtkDataArray creation and insertion
/** This method calls the other factory methods based on extension.

    Caller must delete object after use.
**/

meshBase *meshBase::Create(const std::string &fname) {
  if (fname.find(".vt") != std::string::npos ||
      fname.find(".stl") != std::string::npos) {
    auto *vtkmesh = new vtkMesh(fname);
    vtkmesh->setFileName(fname);
    return vtkmesh;
  } else if (fname.find(".msh") != std::string::npos) {
    std::cout << "Detected file in GMSH format" << std::endl;
    std::cout << "Exporting to VTK ...." << std::endl;
    return exportGmshToVtk(fname);
  } else if (fname.find(".vol") != std::string::npos) {
    std::cout << "Detected file in Netgen .vol format" << std::endl;
    std::cout << "Exporting to VTK ...." << std::endl;
    return exportVolToVtk(fname);
  } else if (fname.find(".pntmesh") != std::string::npos) {
    std::cout << "Detected file in PNTmesh format" << std::endl;
    std::cout << "Processing the file ...." << std::endl;
    return exportPntToVtk(fname);
  } else if (fname.find(".g") != std::string::npos ||
             fname.find(".exo") != std::string::npos) {
    std::cout << "Detected file in Exodus II format" << std::endl;
    std::cout << "Processing the file ...." << std::endl;
    return exportExoToVtk(fname);
  } else {
    std::cout << "mesh files with extension "
              << fname.substr(fname.find_last_of('.')) << " are not supported!"
              << std::endl;
    exit(1);
  }
}

/** Caller must delete object after use.
 **/
meshBase *meshBase::Create(vtkSmartPointer<vtkDataSet> other,
                           const std::string &newname) {
  return new vtkMesh(other, newname);
}

/** Use of this is only valid when mesh has one cell type.
    Caller must delete object after use.
**/
meshBase *meshBase::Create(const std::vector<double> &xCrds,
                           const std::vector<double> &yCrds,
                           const std::vector<double> &zCrds,
                           const std::vector<nemId_t> &elmConn,
                           const int cellType, const std::string &newname) {
  return new vtkMesh(xCrds, yCrds, zCrds, elmConn, cellType, newname);
}

/** Memory is managed by shared pointer, so do not call delete after use.
 **/
std::unique_ptr<meshBase> meshBase::CreateUnique(
    const std::vector<double> &xCrds, const std::vector<double> &yCrds,
    const std::vector<double> &zCrds, const std::vector<nemId_t> &elmConn,
    const int cellType, const std::string &newname) {
  return std::unique_ptr<meshBase>(
      meshBase::Create(xCrds, yCrds, zCrds, elmConn, cellType, newname));
}

/** Memory is managed by shared pointer, so do not call delete after use.
 **/
std::unique_ptr<meshBase> meshBase::CreateUnique(
    vtkSmartPointer<vtkDataSet> other, const std::string &newname) {
  return std::unique_ptr<meshBase>(meshBase::Create(other, newname));
}

/** Memory is managed by shared pointer, so do not call delete after use.
 **/
std::unique_ptr<meshBase> meshBase::CreateUnique(meshBase *mesh) {
  return std::unique_ptr<meshBase>(mesh);
}

/** Memory is managed by shared pointer, so do not call delete after use.
    (be careful with this one!)
**/
std::shared_ptr<meshBase> meshBase::CreateShared(meshBase *mesh) {
  std::shared_ptr<meshBase> sharedMesh;
  sharedMesh.reset(mesh);
  return sharedMesh;
}

/** Memory is managed by shared pointer, so do not call delete after use.
 **/
std::shared_ptr<meshBase> meshBase::CreateShared(
    vtkSmartPointer<vtkDataSet> other, const std::string &newname) {
  std::shared_ptr<meshBase> mesh;
  mesh.reset(meshBase::Create(other, newname));
  return mesh;
}

/** Memory is managed by shared pointer, so do not call delete after use.
 **/
std::shared_ptr<meshBase> meshBase::CreateShared(const std::string &fname) {
  std::shared_ptr<meshBase> mesh;
  mesh.reset(meshBase::Create(fname));
  return mesh;
}

/** Memory is managed by shared pointer, so do not call delete after use.
 **/
std::shared_ptr<meshBase> meshBase::CreateShared(
    const std::vector<double> &xCrds, const std::vector<double> &yCrds,
    const std::vector<double> &zCrds, const std::vector<nemId_t> &elmConn,
    const int cellType, const std::string &newname) {
  std::shared_ptr<meshBase> mesh;
  mesh.reset(meshBase::Create(xCrds, yCrds, zCrds, elmConn, cellType, newname));
  return mesh;
}

/** Memory is managed by shared pointer, so do not call delete after use.
 **/
std::unique_ptr<meshBase> meshBase::CreateUnique(const std::string &fname) {
  return std::unique_ptr<meshBase>(meshBase::Create(fname));
}

/** Caller must delete object after use.
 **/
meshBase *meshBase::generateMesh(const std::string &fname,
                                 const std::string &meshEngine,
                                 meshingParams *params) {
  // TODO: this method should change to incorporate STEP files
  // if ((fname.find(".stl") == -1) && (fname.find(".fms") == -1)) {
  //   std::cerr << "Only CAD files in STL or FMS format are supported"
  //             << std::endl;
  //   exit(1);
  // }

  meshGen *generator = meshGen::Create(fname, meshEngine, params);
  if (generator) {
    int status = generator->createMeshFromSTL(&fname[0u]);
    meshBase *ret = nullptr;
    if (!status) {
      if (meshEngine == "netgen") {
        std::string newname = nemAux::trim_fname(fname, ".vol");
        ret = exportVolToVtk(newname);
      } else if (meshEngine == "gmsh") {
        std::string newname = nemAux::trim_fname(fname, ".msh");
        ret = exportGmshToVtk(newname);
      } else if (meshEngine == "simmetrix") {
        std::string newname = nemAux::trim_fname(fname, ".vtu");
        ret = Create(generator->getDataSet(), newname);
      } else if (meshEngine == "cfmesh") {
        std::string newname = nemAux::trim_fname(fname, ".vtu");
        ret = Create(generator->getDataSet(), newname);
      } else if (meshEngine == "snappyHexMesh") {
        std::string newname = nemAux::trim_fname(fname, ".vtu");
        ret = Create(generator->getDataSet(), newname);
      } else if (meshEngine == "blockMesh") {
        std::string newname = nemAux::trim_fname(fname, ".vtu");
        ret = Create(generator->getDataSet(), newname);
      }
    }
    delete generator;

    if (ret) {
      return ret;
    } else {
      std::cerr << "Mesh Engine " << meshEngine << " not recognized"
                << std::endl;
      exit(1);
    }
  } else {
    std::cerr << "Could not create mesh generator" << std::endl;
    exit(1);
  }
}

/** Caller must delete object after use.
 **/
meshBase *meshBase::stitchMB(const std::vector<meshBase *> &mbObjs) {
  if (!mbObjs.empty()) {
    vtkSmartPointer<vtkAppendFilter> appender =
        vtkSmartPointer<vtkAppendFilter>::New();
    appender->MergePointsOn();
    for (int i = 0; i < mbObjs.size(); ++i) {
      appender->AddInputData(mbObjs[i]->getDataSet());
    }
    appender->Update();
    return meshBase::Create(appender->GetOutput(), "stitched.vtu");
  } else {
    std::cerr << "Nothing to stitch!" << std::endl;
    exit(1);
  }
}

/** This is the shared pointer version of stitchMB. Memory is managed by shared
    pointer, so do not call delete after use.
**/
std::shared_ptr<meshBase> meshBase::stitchMB(
    const std::vector<std::shared_ptr<meshBase>> &_mbObjs) {
  std::vector<meshBase *> mbObjs(_mbObjs.size());
  for (int i = 0; i < mbObjs.size(); ++i) {
    mbObjs[i] = _mbObjs[i].get();
  }
  return meshBase::CreateShared(meshBase::stitchMB(mbObjs));
}

/**
 **/
meshBase *meshBase::extractSelectedCells(meshBase *mesh,
                                         const std::vector<nemId_t> &cellIds) {
  vtkSmartPointer<vtkIdTypeArray> selectionIds =
      vtkSmartPointer<vtkIdTypeArray>::New();
  selectionIds->SetNumberOfComponents(1);
  for (const auto &cellId : cellIds) {
    selectionIds->InsertNextValue(cellId);
  }
  return extractSelectedCells(mesh->getDataSet(), selectionIds);
}

/**
 **/
meshBase *meshBase::extractSelectedCells(
    vtkSmartPointer<vtkDataSet> mesh, vtkSmartPointer<vtkIdTypeArray> cellIds) {
  vtkSmartPointer<vtkSelectionNode> selectionNode =
      vtkSmartPointer<vtkSelectionNode>::New();
  selectionNode->SetFieldType(vtkSelectionNode::CELL);
  selectionNode->SetContentType(vtkSelectionNode::INDICES);
  selectionNode->SetSelectionList(cellIds);
  // create the selection
  vtkSmartPointer<vtkSelection> selection =
      vtkSmartPointer<vtkSelection>::New();
  selection->AddNode(selectionNode);
  // creating extractor
  vtkSmartPointer<vtkExtractSelection> extractSelection =
      vtkSmartPointer<vtkExtractSelection>::New();
  // set stitch surf as input on first port
  extractSelection->SetInputData(0, mesh);
  // set selectionNode as input on second port
  extractSelection->SetInputData(1, selection);
  extractSelection->Update();
  vtkSmartPointer<vtkDataSet> extractedCellMesh =
      vtkDataSet::SafeDownCast(extractSelection->GetOutput());
  meshBase *selectedCellMesh =
      meshBase::Create(extractedCellMesh, "extracted.vtu");
  return selectedCellMesh;
}

/** check for named array in vtk
 **/
int meshBase::IsArrayName(const std::string &name,
                          const bool pointOrCell) const {
  if (!pointOrCell) {
    vtkPointData *pd = dataSet->GetPointData();
    if (pd->GetNumberOfArrays()) {
      for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
        if (name == (pd->GetArrayName(i) ? pd->GetArrayName(i) : "NULL")) {
          return i;
        }
      }
    }
  } else {
    vtkCellData *cd = dataSet->GetCellData();
    if (cd->GetNumberOfArrays()) {
      for (int i = 0; i < cd->GetNumberOfArrays(); ++i) {
        if (name == (cd->GetArrayName(i) ? cd->GetArrayName(i) : "Null")) {
          return i;
        }
      }
    }
  }
  return -1;
}

/** partition mesh into numPartition pieces (static fcn)

    Memory is managed by shared pointer, so do not call delete after use.
**/
std::vector<std::shared_ptr<meshBase>> meshBase::partition(
    const meshBase *mbObj, const int numPartitions) {
  // construct partitioner with meshBase object
  auto *mPart = new meshPartitioner(mbObj);
  if (mPart->partition(numPartitions)) {
    exit(1);
  }
  // initialize vector of meshBase partitions
  std::vector<std::shared_ptr<meshBase>> mbParts(numPartitions);
  for (int i = 0; i < numPartitions; ++i) {
    // define coordinates
    std::vector<std::vector<double>> comp_crds(mbObj->getVertCrds());
    // get partition connectivity and zero index it
    std::vector<int> vtkConn_int(mPart->getConns(i));
    std::vector<nemId_t> vtkConn(vtkConn_int.begin(), vtkConn_int.end());
    for (auto &&it : vtkConn) it -= 1;

    std::string basename(nemAux::trim_fname(mbObj->getFileName(), ""));
    basename += std::to_string(i);
    basename += ".vtu";
    // construct meshBase partition from coordinates and connectivities
    // from partitioner
    mbParts[i] = meshBase::CreateShared(
        mPart->getCrds(i, comp_crds[0]), mPart->getCrds(i, comp_crds[1]),
        mPart->getCrds(i, comp_crds[2]), vtkConn, VTK_TETRA, basename);
    // add partition id to node and cell data of mbPart
    vtkSmartPointer<vtkIntArray> nodePartitionIds =
        vtkSmartPointer<vtkIntArray>::New();
    nodePartitionIds->SetName("NodePartitionIds");
    nodePartitionIds->SetNumberOfComponents(1);
    nodePartitionIds->SetNumberOfTuples(mbParts[i]->getNumberOfPoints());
    nodePartitionIds->FillComponent(0, i);
    mbParts[i]->getDataSet()->GetPointData()->AddArray(nodePartitionIds);

    vtkSmartPointer<vtkIntArray> cellPartitionIds =
        vtkSmartPointer<vtkIntArray>::New();
    cellPartitionIds->SetName("CellPartitionIds");
    cellPartitionIds->SetNumberOfComponents(1);
    cellPartitionIds->SetNumberOfTuples(mbParts[i]->getNumberOfCells());
    cellPartitionIds->FillComponent(0, i);
    mbParts[i]->getDataSet()->GetCellData()->AddArray(cellPartitionIds);

    // add global node index array to partition
    std::map<int, int> partToGlobNodeMap(mPart->getPartToGlobNodeMap(i));
    vtkSmartPointer<vtkIdTypeArray> globalNodeIds =
        vtkSmartPointer<vtkIdTypeArray>::New();
    globalNodeIds->SetName("GlobalNodeIds");
    globalNodeIds->SetNumberOfComponents(1);
    globalNodeIds->SetNumberOfTuples(mbParts[i]->getNumberOfPoints());
    globalNodeIds->SetNumberOfValues(mbParts[i]->getNumberOfPoints());
    auto it = partToGlobNodeMap.begin();
    int idx = 0;
    while (it != partToGlobNodeMap.end()) {
      int globidx = it->second - 1;
      int locidx = it->first - 1;
      globalNodeIds->SetTuple1(idx, globidx);
      mbParts[i]->globToPartNodeMap[globidx] = locidx;
      mbParts[i]->partToGlobNodeMap[locidx] = globidx;
      ++idx;
      ++it;
    }
    mbParts[i]->getDataSet()->GetPointData()->AddArray(globalNodeIds);
    // add global cell index array to partition
    std::map<int, int> partToGlobCellMap(mPart->getPartToGlobElmMap(i));
    // vtkSmartPointer<vtkIdTypeArray> globalCellIds =
    // vtkSmartPointer<vtkIdTypeArray>::New();
    // globalCellIds->SetName("GlobalCellIds");
    // globalCellIds->SetNumberOfComponents(1);
    // globalCellIds->SetNumberOfTuples(mbPart->getNumberOfCells());
    // globalCellIds->SetNumberOfValues(mbPart->getNumberOfCells());
    it = partToGlobCellMap.begin();
    idx = 0;
    while (it != partToGlobCellMap.end()) {
      int globidx = it->second - 1;
      int locidx = it->first - 1;
      //    globalCellIds->SetTuple1(idx,globidx);
      mbParts[i]->globToPartCellMap[globidx] = locidx;
      mbParts[i]->partToGlobCellMap[locidx] = globidx;
      ++idx;
      ++it;
    }
    // mbPart->getDataSet()->GetCellData()->AddArray(globalCellIds);
    mbParts[i]->write();
    // mbParts[i] = mbPart;
  }
  delete mPart;
  mPart = nullptr;
  return mbParts;
}

/** integrated data is available per cell as vtk cell data after operation
 **/
std::vector<std::vector<double>> meshBase::integrateOverMesh(
    const std::vector<int> &arrayIDs) {
  std::unique_ptr<GaussCubature> cubature =
      GaussCubature::CreateUnique(dataSet, arrayIDs);
  return cubature->integrateOverAllCells();
}

/**
 **/
void meshBase::generateSizeField(const std::string &method, int arrayID,
                                 double dev_mult, bool maxIsmin,
                                 double sizeFactor, int order) {
  std::cout << "Size Factor = " << sizeFactor << std::endl;
  std::unique_ptr<NEM::ADP::SizeFieldBase> sfobj =
      NEM::ADP::SizeFieldBase::CreateUnique(dataSet, method, arrayID, dev_mult,
                                            maxIsmin, sizeFactor, order);
  sfobj->setSizeFactor(sizeFactor);
  sfobj->computeSizeField(dataSet->GetPointData()->GetArray(arrayID));
}

/**
 **/
meshBase *meshBase::exportGmshToVtk(const std::string &fname) {
  std::ifstream meshStream(fname);
  if (!meshStream.good()) {
    std::cout << "Error opening file " << fname << std::endl;
    exit(1);
  }

  bool warning = true;

  std::string line;
  int numPoints = 0, numCells = 0, numPhysGrps = 0;
  bool fndPhyGrp = false;
  std::vector<int> physGrpDims;
  std::map<int, std::string> physGrpIdName;
  std::vector<std::vector<std::vector<double>>> pointData;
  std::vector<std::vector<std::vector<double>>> cellData;
  std::vector<std::vector<double>> cellPhysGrpIds;
  std::vector<std::string> pointDataNames;
  std::vector<std::string> cellDataNames;
  std::vector<std::string> fieldDataNames;
  std::vector<std::string> fieldData;

  // declare points to be pushed into dataSet_tmp
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToDouble();
  // declare dataSet_tmp which will be associated to output vtkMesh
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  // map to hold true index of points (gmsh allows non-contiguous ordering)
  std::map<int, int> trueIndex;

  while (getline(meshStream, line)) {
    if (line.find("$PhysicalNames") != -1) {
      fndPhyGrp = true;
      getline(meshStream, line);
      std::stringstream ss(line);
      ss >> numPhysGrps;
      std::cout << "Found " << numPhysGrps << " physical groups!\n";
      for (int i = 0; i < numPhysGrps; ++i) {
        int grpDim, grpId;
        std::string grpName;
        getline(meshStream, line);
        std::stringstream ss(line);
        ss >> grpDim >> grpId >> grpName;
        grpName.erase(std::remove(grpName.begin(), grpName.end(), '\"'),
                      grpName.end());
        /*std::cout << "Group Dim = " << grpDim
                    << "\nGroup Id = " << grpId
                    << "\nGroup Name = " << grpName
                    << std::endl;*/
        physGrpDims.push_back(grpDim);
        physGrpIdName[grpId] = grpName;
      }
    }

    if (line.find("$Nodes") != -1) {
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
        ss >> id >> std::setprecision(16) >> x >> y >> z;
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

    if (line.find("$Elements") != -1) {
      getline(meshStream, line);
      // std::cout << "line = " << line << std::endl;
      std::stringstream ss(line);
      ss >> numCells;
      int id, type, numTags;
      // double tmp2[1];
      // allocate space for cell connectivities
      dataSet_tmp->Allocate(numCells);
      for (int i = 0; i < numCells; ++i) {
        getline(meshStream, line);
        std::stringstream ss(line);
        ss >> id >> type >> numTags;
        vtkSmartPointer<vtkIdList> vtkCellIds =
            vtkSmartPointer<vtkIdList>::New();
        if (type == 2) {
          int tmp;
          if (!fndPhyGrp) {
            for (int j = 0; j < numTags; ++j) ss >> tmp;
          } else {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags - 1; ++j) ss >> tmp;
          }
          for (int j = 0; j < 3; ++j) {
            ss >> tmp;
            // insert connectivities for cell into cellIds container
            vtkCellIds->InsertNextId(trueIndex[tmp]);  //-1);
          }
          // insert connectivities for triangle elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_TRIANGLE, vtkCellIds);
        } else if (type == 4) {
          int tmp;
          if (!fndPhyGrp) {
            for (int j = 0; j < numTags; ++j) ss >> tmp;
          } else {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags - 1; ++j) ss >> tmp;
          }
          for (int j = 0; j < 4; ++j) {
            ss >> tmp;
            // insert connectivities for cell into cellIds container
            vtkCellIds->InsertNextId(trueIndex[tmp]);  //-1);
          }
          // insert connectivities for tet elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_TETRA, vtkCellIds);
        } else if (type == 3) {
          int tmp;
          if (!fndPhyGrp) {
            for (int j = 0; j < numTags; ++j) ss >> tmp;
          } else {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags - 1; ++j) ss >> tmp;
          }
          for (int j = 0; j < 4; ++j) {
            ss >> tmp;
            // insert connectivities for cell into cellIds container
            vtkCellIds->InsertNextId(trueIndex[tmp]);  //-1);
          }
          // insert connectivities for tet elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_QUAD, vtkCellIds);
        } else if (type == 5) {
          int tmp;
          if (!fndPhyGrp) {
            for (int j = 0; j < numTags; ++j) ss >> tmp;
          } else {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags - 1; ++j) ss >> tmp;
          }
          for (int j = 0; j < 8; ++j) {
            ss >> tmp;
            // insert connectivities for cell into cellIds container
            vtkCellIds->InsertNextId(trueIndex[tmp]);  //-1);
          }
          // insert connectivities for hex elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_HEXAHEDRON, vtkCellIds);
        } else if (type == 6) {
          int tmp;
          if (!fndPhyGrp) {
            for (int j = 0; j < numTags; ++j) ss >> tmp;
          } else {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags - 1; ++j) ss >> tmp;
          }
          for (int j = 0; j < 6; ++j) {
            ss >> tmp;
            // insert connectivities for cell into cellIds container
            vtkCellIds->InsertNextId(trueIndex[tmp]);  //-1);
          }
          // insert connectivities for wedge/prism elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_WEDGE, vtkCellIds);
        } else if (type == 11) {
          int tmp;
          if (!fndPhyGrp) {
            for (int j = 0; j < numTags; ++j) ss >> tmp;
          } else {
            std::vector<double> physGrpId(1);
            ss >> physGrpId[0];
            cellPhysGrpIds.push_back(physGrpId);
            for (int j = 0; j < numTags - 1; ++j) ss >> tmp;
          }
          for (int j = 0; j < 10; ++j) {
            ss >> tmp;
            // insert connectivities for cell into cellIds container
            vtkCellIds->InsertNextId(trueIndex[tmp]);  //-1);
          }
          // insert connectivities for tet elements into dataSet
          dataSet_tmp->InsertNextCell(VTK_QUADRATIC_TETRA, vtkCellIds);
        } else {
          if (warning) {
            std::cout
                << "Warning: Only triangular, quadrilateral, "
                   "tetrahedral, hexahedral, wedge, and quadratic tetrahedral"
                   "elements are supported, "
                << "everything else is ignored! " << std::endl;
            warning = false;
            // exit(1);
          }
        }
      }
    }

    if (line.find("$NodeData") != -1) {
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
      dataname.erase(std::remove(dataname.begin(), dataname.end(), '\"'),
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
      int dt, dim, numFields, tmp;
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

    if (line.find("$ElementData") != -1) {
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
      dataname.erase(std::remove(dataname.begin(), dataname.end(), '\"'),
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
      int dt, dim, numFields, tmp;
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

  if (fndPhyGrp) {
    cellDataNames.push_back("PhysGrpId");
    cellData.push_back(cellPhysGrpIds);
  }

  vtkMesh *vtkmesh = new vtkMesh();
  // vtkmesh->dataSet = dataSet_tmp->NewInstance();
  // vtkmesh->dataSet->DeepCopy(dataSet_tmp);
  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  if (numCells == 0) std::cerr << "Warning: No cells were found in the mesh!\n";
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();
  if (numPoints == 0)
    std::cerr << "Warning: No points were found in the mesh!\n";

  if (numPoints > 0) {
    for (int i = 0; i < pointData.size(); ++i)
      vtkmesh->setPointDataArray(&(pointDataNames[i])[0u], pointData[i]);
  }

  if (numCells > 0) {
    for (int i = 0; i < cellData.size(); ++i)
      vtkmesh->setCellDataArray(&(cellDataNames[i])[0u], cellData[i]);
  }

  // vtkmesh->setFileName(nemAux::trim_fname(fname, ".vtu"));
  // vtkmesh->write();
  std::cout << "vtkMesh constructed" << std::endl;

  return vtkmesh;
}

/**
 **/
meshBase *meshBase::exportVolToVtk(const std::string &fname) {
#ifdef HAVE_NGEN
  nglib::Ng_Mesh *Ngmesh;
  nglib::Ng_Init();
  Ngmesh = nglib::Ng_NewMesh();

  int status = nglib::Ng_MergeMesh(Ngmesh, &fname[0u]);
  if (status) {
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
  for (int i = 1; i <= numNgPoints; ++i) {
    double point[3];
    nglib::Ng_GetPoint(Ngmesh, i, point);
    // insert point i
    points->SetPoint(i - 1, point);
  }
  // inserting point array into dataSet_tmp
  dataSet_tmp->SetPoints(points);

  // allocating space for cell connectivities
  dataSet_tmp->Allocate(numSurfCells + numVolCells);

  for (int i = 1; i <= numSurfCells; ++i) {
    vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
    int tri[3];
    nglib::Ng_GetSurfaceElement(Ngmesh, i, tri);
    for (int j = 0; j < 3; ++j) {
      // insert connectivies for cell into cellIds container
      vtkcellIds->InsertNextId(tri[j] - 1);
    }
    // insert connectivies for triangle elements into dataSet
    dataSet_tmp->InsertNextCell(VTK_TRIANGLE, vtkcellIds);
  }

  for (int i = 1; i <= numVolCells; ++i) {
    vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
    int tet[4];
    nglib::Ng_GetVolumeElement(Ngmesh, i, tet);
    for (int j = 0; j < 4; ++j) {
      // insert connectivies for cell into cellIds container
      vtkcellIds->InsertNextId(tet[j] - 1);
    }
    // insert connectivies for triangle elements into dataSet
    dataSet_tmp->InsertNextCell(VTK_TETRA, vtkcellIds);
  }

  vtkMesh *vtkmesh = new vtkMesh();
  // vtkmesh->dataSet = dataSet_tmp->NewInstance();//
  // vtkmesh->dataSet->DeepCopy(dataSet_tmp);//vtkDataSet::SafeDownCast(dataSet_tmp));
  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();

  vtkmesh->setFileName(nemAux::trim_fname(fname, ".vtu"));
  std::cout << "vtkMesh constructed" << std::endl;

  if (Ngmesh) nglib::Ng_DeleteMesh(Ngmesh);
  nglib::Ng_Exit();
  return vtkmesh;
#else
  std::cerr << "ENABLE_NETGEN is not used during build."
            << " Build with ENABLE_NETGEN to use this function." << std::endl;
  exit(1);
#endif
}

/** exports pntMesh to vtk format
 **/
meshBase *meshBase::exportPntToVtk(const std::string &fname) {
  PNTMesh::pntMesh *pMesh;
  pMesh = new PNTMesh::pntMesh(fname);

  vtkMesh *vtkmesh = new vtkMesh();

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
  int numVolCells = pMesh->getNumberOfCells();

  // allocate memory for points
  points->SetNumberOfPoints(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    std::vector<double> point;
    point = pMesh->getPointCrd(i);
    // insert point i
    points->SetPoint(i, &point[0]);
  }
  // inserting point array into dataSet_tmp
  dataSet_tmp->SetPoints(points);

  // allocating space for cell connectivities
  dataSet_tmp->Allocate(numVolCells);
  for (int i = 0; i < numVolCells; ++i) {
    std::cout << "Element " << i << std::endl;
    vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
    std::vector<int> conn;
    VTKCellType vct =
        pMesh->getVtkCellTag(pMesh->getElmType(i), pMesh->getElmOrder(i));
    conn = pMesh->getElmConn(i, vct);
    for (int j = 0; j < conn.size(); ++j) {
      // insert connectivies for cell into cellIds container
      vtkcellIds->InsertNextId(conn[j]);
    }
    // insert connectivies
    dataSet_tmp->InsertNextCell(vct, vtkcellIds);
  }

  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();

  //
  std::cout << "Trimmed name = " << nemAux::trim_fname(fname, ".vtu")
            << std::endl;
  vtkmesh->setFileName(nemAux::trim_fname(fname, ".vtu"));
  // vtkmesh->write();
  std::cout << "vtkMesh constructed" << std::endl;

  return vtkmesh;
}

/** exports exodusII to vtk format
 **/
meshBase *meshBase::exportExoToVtk(const std::string &fname) {
  // opening the file
  int CPU_word_size, IO_word_size;
  int fid, _exErr;
  float version;
  CPU_word_size = sizeof(float);
  IO_word_size = 0;
  _exErr = 0;

  /* open EXODUS II files */
  fid =
      ex_open(fname.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem opening file " + fname + "\n");

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
  char cdum;
  _exErr = ex_inquire(fid, EX_INQ_API_VERS, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II API version is " << fdum << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_DB_VERS, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Exodus II Database version is " << fdum << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_DIM, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  std::cout << "Number of coordinate dimensions is " << num_props << std::endl;
  if (num_props != 3)
    NEM::MSH::EXOMesh::wrnErrMsg(-1, "Only 3D mesh data is supported!\n");

  _exErr = ex_inquire(fid, EX_INQ_NODES, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numPoints = num_props;
  std::cout << "Number of points " << numPoints << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_ELEM, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numVolCells = num_props;
  std::cout << "Number of elements " << numVolCells << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_ELEM_BLK, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numElmBlk = num_props;
  std::cout << "Number of element blocks " << numElmBlk << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_NODE_SETS, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numNdeSet = num_props;
  std::cout << "Number of node sets " << numNdeSet << std::endl;

  _exErr = ex_inquire(fid, EX_INQ_SIDE_SETS, &num_props, &fdum, &cdum);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading file contents.\n");
  numSideSet = num_props;
  std::cout << "Number of side sets " << numSideSet << std::endl;

  // nodal coordinates
  std::vector<float> x, y, z;
  x.resize(numPoints, 0);
  y.resize(numPoints, 0);
  z.resize(numPoints, 0);
  _exErr = ex_get_coord(fid, &x[0], &y[0], &z[0]);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem reading nodal coordinates.\n");

  // allocate memory for points
  points->SetNumberOfPoints(numPoints);
  for (int i = 0; i < numPoints; ++i) {
    std::vector<double> pt = {x[i], y[i], z[i]};
    points->SetPoint(i, &pt[0]);
  }
  // inserting point array into dataSet_tmp
  dataSet_tmp->SetPoints(points);

  // Connectivities
  // allocating space for cell connectivities
  dataSet_tmp->Allocate(numVolCells);
  // read element blocks
  for (int iEB = 1; iEB <= numElmBlk; iEB++) {
    int num_el_in_blk, num_nod_per_el, num_attr /*, *connect*/;
    // float *attrib;
    char elem_type[MAX_STR_LENGTH + 1];
    // read element block parameters
    _exErr = ex_get_elem_block(fid, iEB, elem_type, &num_el_in_blk,
                               &num_nod_per_el, &num_attr);
    NEM::MSH::EXOMesh::wrnErrMsg(_exErr,
                                 "Problem reading element block parameters.\n");
    // read element connectivity
    std::vector<int> conn;
    conn.resize(num_el_in_blk * num_nod_per_el, 0);
    _exErr = ex_get_elem_conn(fid, iEB, &conn[0]);
    NEM::MSH::EXOMesh::wrnErrMsg(
        _exErr, "Problem reading element block connectivities.\n");
    // read element block attributes
    // std::vector<float> attr;
    // attr.resize(0.,num_el_in_blk*num_nod_per_el);
    //_exErr = ex_get_elem_attr (fid, iEB, &attrib[0]);
    // EXOMesh::wrnErrMsg(_exErr, "Problem reading element block
    // attributes.\n");
    for (int iEl = 0; iEl < num_el_in_blk; ++iEl) {
      vtkSmartPointer<vtkIdList> vtkcellIds = vtkSmartPointer<vtkIdList>::New();
      VTKCellType vct =
          NEM::MSH::EXOMesh::e2vEMap(NEM::MSH::EXOMesh::elmTypeNum(elem_type));
      for (int jc = iEl * num_nod_per_el; jc < (iEl + 1) * num_nod_per_el;
           jc++) {
        // insert connectivities for cell into cellIds container
        vtkcellIds->InsertNextId(conn[jc] - 1);
      }
      // insert connectivities
      dataSet_tmp->InsertNextCell(vct, vtkcellIds);
    }
  }

  std::cout << "Trimmed name = " << nemAux::trim_fname(fname, ".vtu")
            << std::endl;

  vtkMesh *vtkmesh = new vtkMesh();
  vtkmesh->dataSet = dataSet_tmp;
  vtkmesh->numCells = vtkmesh->dataSet->GetNumberOfCells();
  vtkmesh->numPoints = vtkmesh->dataSet->GetNumberOfPoints();
  vtkmesh->setFileName(nemAux::trim_fname(fname, ".vtu"));
  // vtkmesh->write();
  std::cout << "vtkMesh constructed" << std::endl;

  // closing the file
  _exErr = ex_close(fid);
  NEM::MSH::EXOMesh::wrnErrMsg(_exErr, "Problem closing the exodusII file.");

  return vtkmesh;
}

/** convert to gmsh format without data
 **/
void meshBase::writeMSH(std::ofstream &outputStream) {
  if (!outputStream.good()) {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }

  if (!dataSet) {
    std::cout << "No data to write" << std::endl;
    exit(1);
  }
  // ---------  writing gmsh header ----------- //
  outputStream << "$MeshFormat" << std::endl
               << "2.2 0 8" << std::endl
               << "$EndMeshFormat" << std::endl;

  // -------- ensure all cell types are tri/tet or below -------------- //
  int type_id;
  for (int i = 0; i < numCells; i++) {
    type_id = dataSet->GetCellType(i);
    if (!(type_id == 3 || type_id == 5 || type_id == 10 || type_id == 9)) {
      std::cout << "Error: Only tetrahedral, triangular, and quadrilateral"
                << " meshes can be written to gmsh format" << std::endl;
      exit(3);
    }
  }

  // ------------------------ write point coords -------------------------- //
  outputStream << "$Nodes" << std::endl << numPoints << std::endl;
  for (int i = 0; i < numPoints; ++i) {
    std::vector<double> pntcrds = getPoint(i);
    outputStream << i + 1 << " " << pntcrds[0] << " " << pntcrds[1] << " "
                 << pntcrds[2] << " " << std::endl;
  }
  outputStream << "$EndNodes" << std::endl;

  // ------------- write element type and connectivity --------------------- //
  outputStream << "$Elements" << std::endl << numCells << std::endl;
  for (int i = 0; i < numCells; ++i) {
    vtkIdList *point_ids = dataSet->GetCell(i)->GetPointIds();
    int numComponent = point_ids->GetNumberOfIds();
    type_id = dataSet->GetCell(i)->GetCellType();
    outputStream << i + 1 << " ";
    switch (numComponent) {
      case 2: {
        outputStream << 1 << " " << 2 << " " << 1 << " " << 1 << " ";
        break;
      }
      case 3: {
        outputStream << 2 << " " << 2 << " " << 1 << " " << 1 << " ";
        break;
      }
      case 4: {
        if (type_id == 10)  // tetra
        {
          outputStream << 4 << " " << 2 << " " << 1 << " " << 1 << " ";
          break;
        } else if (type_id == 9)  // quad
        {
          outputStream << 3 << " " << 2 << " " << 1 << " " << 1 << " ";
          break;
        }
      }

      default: {
        std::cerr << "Components in cell should be <= 4" << std::endl;
        exit(1);
      }
    }
    for (int j = 0; j < numComponent; ++j)
      outputStream << point_ids->GetId(j) + 1 << " ";
    outputStream << std::endl;
  }
  outputStream << "$EndElements" << std::endl;
}

/** convert to gmsh format with specified point or cell data
 **/
void meshBase::writeMSH(std::ofstream &outputStream,
                        const std::string &pointOrCell, int arrayID) {
  // write points and cells
  writeMSH(outputStream);

  if (!outputStream.good()) {
    std::cout << "Output file stream is bad" << std::endl;
    exit(1);
  }

  if (!pointOrCell.compare("point")) {
    // ---------------------------- write point data -------------------------
    // //
    vtkPointData *pointData = dataSet->GetPointData();
    if (pointData) {
      int numArr = pointData->GetNumberOfArrays();
      if (arrayID >= numArr) {
        std::cout << "ERROR: arrayID is out of bounds" << std::endl;
        std::cout << "There are " << numArr << " point data arrays"
                  << std::endl;
        exit(1);
      } else if (numArr < 1) {
        std::cout << "no point data found" << std::endl;
        exit(1);
      }
      vtkDataArray *da = pointData->GetArray(arrayID);
      if (da) {
        int numComponent = da->GetNumberOfComponents();
        int numTuple = da->GetNumberOfTuples();
        std::string tmpname = "PointArray";
        tmpname += std::to_string(arrayID);
        outputStream << "$NodeData" << std::endl
                     << 1 << std::endl  // 1 string tag
                     << "\""
                     << (pointData->GetArrayName(arrayID)
                             ? pointData->GetArrayName(arrayID)
                             : tmpname)  // name of view
                     << "\"" << std::endl
                     << 0 << std::endl  // 0 real tag
                     << 3 << std::endl  // 3 int tags (dt index, dim of field,
                                        // number of fields)
                     << 0 << std::endl  // dt index
                     << numComponent << std::endl  // dim of field
                     << numTuple << std::endl;     // number of fields
        for (int j = 0; j < numTuple; ++j) {
          double *data = da->GetTuple(j);
          outputStream << j + 1 << " ";
          for (int k = 0; k < numComponent; ++k) {
            outputStream << data[k] << " ";
          }
          outputStream << std::endl;
        }
        outputStream << "$EndNodeData" << std::endl;
      }
    }
  } else if (!pointOrCell.compare("cell")) {
    // -------------------------- write cell data ----------------------------
    // //

    vtkCellData *cellData = dataSet->GetCellData();
    if (cellData) {
      int numArr = cellData->GetNumberOfArrays();
      if (arrayID >= numArr) {
        std::cout << "ERROR: arrayID is out of bounds" << std::endl;
        std::cout << "There are " << numArr << " cell data arrays" << std::endl;
        exit(1);
      } else if (numArr < 1) {
        std::cout << "no cell data found" << std::endl;
        exit(1);
      }
      vtkDataArray *da = cellData->GetArray(arrayID);
      if (da) {
        int numComponent = da->GetNumberOfComponents();
        int numTuple = da->GetNumberOfTuples();
        std::string tmpname = "CellArray";
        tmpname += std::to_string(arrayID);
        outputStream << "$ElementData" << std::endl
                     << 1 << std::endl  // 1 string tag
                     << "\""
                     << (cellData->GetArrayName(arrayID)
                             ? cellData->GetArrayName(arrayID)
                             : tmpname)  // name of view
                     << "\"" << std::endl
                     << 0 << std::endl  // 0 real tag
                     << 3 << std::endl  // 3 int tags (dt index, dim of field,
                                        // number of fields)
                     << 0 << std::endl  // dt index
                     << numComponent << std::endl  // dim of field
                     << numTuple << std::endl;     // number of fields
        for (int j = 0; j < numTuple; ++j) {
          double *data = da->GetTuple(j);
          outputStream << j + 1 << " ";
          for (int k = 0; k < numComponent; ++k) {
            outputStream << data[k] << " ";
          }
          outputStream << std::endl;
        }
        outputStream << "$EndElementData" << std::endl;
      }
    }
  }
}

/** convert to gmsh format with specified point or cell data for
 **/
void meshBase::writeMSH(std::ofstream &outputStream,
                        const std::string &pointOrCell, int arrayID,
                        bool onlyVol) {
  if (!outputStream.good()) {
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
  for (int i = 0; i < numCells; i++) {
    int type_id = dataSet->GetCellType(i);
    if (!(type_id == 3 || type_id == 5 || type_id == 10 || type_id == 9)) {
      std::cout << "Error: Only tetrahedral, triangular, and quadrilateral"
                << " meshes can be written to gmsh format" << std::endl;
      exit(3);
    }
    if (!(type_id == 10)) num_bad += 1;
  }

  // ------------------------ write point coords -------------------------- //
  outputStream << "$Nodes" << std::endl << numPoints << std::endl;
  for (int i = 0; i < numPoints; ++i) {
    std::vector<double> pntcrds = getPoint(i);
    outputStream << i + 1 << " ";
    outputStream << std::setprecision(16)
                 << pntcrds[0] << " " << pntcrds[1] << " " << pntcrds[2]
                 << " " << std::endl;
  }
  outputStream << "$EndNodes" << std::endl;

  // ------------- write element type and connectivity --------------------- //
  outputStream << "$Elements" << std::endl << numCells - num_bad << std::endl;
  // int k = 0;
  for (int i = 0; i < numCells; ++i) {
    vtkIdList *point_ids = dataSet->GetCell(i)->GetPointIds();
    int numComponent = point_ids->GetNumberOfIds();
    int type_id = dataSet->GetCellType(i);
    if (type_id == 10) {
      outputStream << i + 1 << " ";
      switch (numComponent) {
        case 2: {
          break;
        }
        case 3: {
          outputStream << 2 << " " << 2 << " " << 1 << " " << 1 << " ";
          break;
        }
        case 4: {
          outputStream << 4 << " " << 2 << " " << 1 << " " << 1 << " ";
          break;
        }

        default: {
          std::cerr << "Components in cell should be less than 4" << std::endl;
          exit(1);
        }
      }
      for (int j = 0; j < numComponent; ++j)
        outputStream << point_ids->GetId(j) + 1 << " ";
      outputStream << std::endl;
      // k+=1;
    }
  }
  outputStream << "$EndElements" << std::endl;
  // -------------------------- write cell data ---------------------------- //
  vtkCellData *cellData = dataSet->GetCellData();
  vtkDataArray *da = cellData->GetArray(arrayID);
  if (da) {
    std::string tmpname = "CellArray";
    tmpname += std::to_string(arrayID);
    outputStream
        << "$ElementData" << std::endl
        << 1 << std::endl  // 1 string tag
        << "\""
        << (cellData->GetArrayName(arrayID) ? cellData->GetArrayName(arrayID)
                                            : tmpname)  // name of view
        << "\"" << std::endl
        << 0 << std::endl  // 0 real tag
        << 3
        << std::endl  // 3 int tags (dt index, dim of field, number of fields)
        << 0 << std::endl                    // dt index
        << 1 << std::endl                    // dim of field
        << numCells - num_bad << std::endl;  // number of fields
    for (int j = 0; j < numCells; ++j) {
      int type_id = dataSet->GetCellType(j);
      if (type_id == 10) {
        double *data = da->GetTuple(j);
        outputStream << j + 1 << " ";
        outputStream << data[0] << " ";
        outputStream << std::endl;
      }
    }
    outputStream << "$EndElementData" << std::endl;
  }
}

/**
 **/
void writePatchMap(const std::string &mapFile,
                   const std::map<int, int> &patchMap) {
  std::ofstream outputStream(mapFile);
  if (!outputStream.good()) {
    std::cout << "Error opening file " << mapFile << std::endl;
    exit(1);
  }
  writePatchMap(outputStream, patchMap);
}

/**
 **/
void writePatchMap(std::ofstream &outputStream,
                   const std::map<int, int> &patchMap) {
  outputStream << patchMap.size() << std::endl;
  outputStream << patchMap.size() << std::endl;
  auto it = patchMap.begin();
  int normPatchNo = 1;
  while (it != patchMap.end()) {
    for (int i = 0; i < 2; ++i) {
      outputStream << std::setw(2) << std::left << it->first << " ";
    }
    outputStream << std::setw(2) << std::left << normPatchNo << " ";
    outputStream << std::endl;
    ++it;
    normPatchNo++;
  }
}

/** for rocstar restart hack through rflupart/prep
 **/
void meshBase::writeCobalt(meshBase *surfWithPatches,
                           const std::string &mapFile,
                           std::ofstream &outputStream) {
  if (!surfWithPatches) {
    std::cout << "surface mesh is empty!" << std::endl;
    exit(1);
  }
  if (surfWithPatches->IsArrayName("patchNo", 1) == -1) {
    std::cout << "surface mesh must have patchNo cell array" << std::endl;
    exit(1);
  }
  vtkSmartPointer<vtkIdList> facePtIds;
  vtkSmartPointer<vtkIdList> sharedCellPtIds =
      vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkGenericCell> genCell1 =
      vtkSmartPointer<vtkGenericCell>::New();
  vtkSmartPointer<vtkGenericCell> genCell2 =
      vtkSmartPointer<vtkGenericCell>::New();
  std::map<std::vector<nemId_t>, std::pair<nemId_t, nemId_t>,
           sortNemId_tVec_compare>
      faceMap;
  // building cell locator for looking up patch number in remeshed surface mesh
  vtkSmartPointer<vtkStaticCellLocator> surfCellLocator =
    surfWithPatches->buildStaticCellLocator();
  // maximum number of vertices per face (to be found in proceeding loop)
  int nVerticesPerFaceMax = 0;
  // maximum number of faces per cell (to be found in proceeding loop)
  int nFacesPerCellMax = 0;

  for (nemId_t i = 0; i < this->getNumberOfCells(); ++i) {
    // get cell i
    dataSet->GetCell(i, genCell1);
    // get faces, find cells sharing it. if no cell shares it,
    // use the locator of the surfWithPatches to find the patch number
    int numFaces = genCell1->GetNumberOfFaces();
    nFacesPerCellMax =
        (nFacesPerCellMax < numFaces ? numFaces : nFacesPerCellMax);
    for (int j = 0; j < numFaces; ++j) {
      vtkCell *face = genCell1->GetFace(j);
      // bool shared = false;
      vtkIdType numVerts = face->GetNumberOfPoints();
      nVerticesPerFaceMax =
          (nVerticesPerFaceMax < numVerts ? numVerts : nVerticesPerFaceMax);
      facePtIds = face->GetPointIds();
      dataSet->GetCellNeighbors(i, facePtIds, sharedCellPtIds);
      std::vector<nemId_t> facePntIds(numVerts);
      for (vtkIdType k = 0; k < numVerts; ++k) {
        facePntIds[k] = face->GetPointId(k) + 1;
      }
      // std::cout << "sharedCellPtIds->GetNumberOfIds() = " <<
      // sharedCellPtIds->GetNumberOfIds() << std::endl;
      if (sharedCellPtIds->GetNumberOfIds()) {
        faceMap.insert(
            std::pair<std::vector<nemId_t>, std::pair<nemId_t, nemId_t>>(
                facePntIds,
                std::make_pair(i + 1, sharedCellPtIds->GetId(0) + 1)));
      } else {
        double p1[3], p2[3], p3[3];
        face->GetPoints()->GetPoint(0, p1);
        face->GetPoints()->GetPoint(1, p2);
        face->GetPoints()->GetPoint(2, p3);
        double faceCenter[3];
        for (vtkIdType k = 0; k < numVerts; ++k) {
          faceCenter[k] = (p1[k] + p2[k] + p3[k]) / 3.0;
        }
        vtkIdType closestCellId;
        int subId;
        double minDist2;
        double closestPoint[3];
        // find closest point and closest cell to faceCenter
        surfCellLocator->FindClosestPoint(faceCenter, closestPoint, genCell2,
                                          closestCellId, subId, minDist2);
        double patchNo[1];
        surfWithPatches->getDataSet()
            ->GetCellData()
            ->GetArray("patchNo")
            ->GetTuple(closestCellId, patchNo);
        faceMap.insert(
            std::pair<std::vector<nemId_t>, std::pair<nemId_t, nemId_t>>(
                facePntIds, std::make_pair(i + 1, -1 * patchNo[0])));
      }
    }
  }

  std::map<int, int> patchMap;
  for (int i = 0; i < surfWithPatches->getNumberOfCells(); ++i) {
    double patchNo[1];
    surfWithPatches->getDataSet()->GetCellData()->GetArray("patchNo")->GetTuple(
        i, patchNo);
    patchMap.insert(std::pair<int, int>(patchNo[0], i));
  }

  // write patch mapping file
  writePatchMap(mapFile, patchMap);
  // write cobalt file
  outputStream << 3 << "   " << 1 << "  " << patchMap.size() << std::endl;
  outputStream << this->getNumberOfPoints() << " " << faceMap.size() << " "
               << this->getNumberOfCells() << " " << nVerticesPerFaceMax << " "
               << nFacesPerCellMax << std::endl;
  for (int i = 0; i < this->getNumberOfPoints(); ++i) {
    std::vector<double> pnt(this->getPoint(i));
    outputStream << std::setw(21) << std::fixed << std::setprecision(15)
                 << pnt[0] << "   " << pnt[1] << "   " << pnt[2] << std::endl;
  }

  auto it = faceMap.begin();
  while (it != faceMap.end()) {
    outputStream << it->first.size() << " ";
    auto faceIdIter = it->first.begin();
    while (faceIdIter != it->first.end()) {
      outputStream << *faceIdIter << " ";
      ++faceIdIter;
    }
    outputStream << it->second.first << " " << it->second.second << std::endl;
    ++it;
  }
}

/** for rocstar restart hack through rflupart/prep
 **/
void meshBase::writeCobalt(meshBase *surfWithPatches,
                           const std::string &mapFile,
                           const std::string &ofname) {
  std::ofstream outputStream(ofname);
  if (!outputStream.good()) {
    std::cout << "Cannot open file " << ofname << std::endl;
    exit(1);
  }
  writeCobalt(surfWithPatches, mapFile, outputStream);
}

/**
 **/
void meshBase::writeMSH(const std::string &fname,
                        const std::string &pointOrCell, int arrayID,
                        bool onlyVol) {
  std::ofstream outputStream(fname.c_str());
  writeMSH(outputStream, pointOrCell, arrayID, onlyVol);
}

/**
 **/
void meshBase::writeMSH(const std::string &fname) {
  std::ofstream outputStream(fname.c_str());
  writeMSH(outputStream);
}

/**
 **/
void meshBase::writeMSH(const std::string &fname,
                        const std::string &pointOrCell, int arrayID) {
  std::ofstream outputStream(fname.c_str());
  writeMSH(outputStream, pointOrCell, arrayID);
}

/** edge_scale is for uniform refinement and is ignored in calls where
    method is "gradient", "value", etc. instead of "uniform"
**/
void meshBase::refineMesh(const std::string &method, int arrayID,
                          double dev_mult, bool maxIsmin, double edge_scale,
                          const std::string &ofname, bool transferData,
                          double sizeFactor, bool constrainBoundary) {
  std::unique_ptr<NEM::ADP::Refine> refineobj =
      std::unique_ptr<NEM::ADP::Refine>(
          new NEM::ADP::Refine(this, method, arrayID, dev_mult, maxIsmin,
                               edge_scale, ofname, sizeFactor));
  refineobj->run(transferData, constrainBoundary);
}

/**
 **/
void meshBase::refineMesh(const std::string &method, int arrayID, int _order,
                          const std::string &ofname, bool transferData) {
  std::unique_ptr<NEM::ADP::Refine> refineobj =
      std::unique_ptr<NEM::ADP::Refine>(new NEM::ADP::Refine(
          this, method, arrayID, 0.0, false, 0.0, ofname, 1.0, _order));
  refineobj->run(transferData);
}

/** edge_scale is for uniform refinement and is ignored in calls where
    method is "gradient", "value", etc. instead of "uniform"
**/
void meshBase::refineMesh(const std::string &method,
                          const std::string &arrayName, double dev_mult,
                          bool maxIsmin, double edge_scale,
                          const std::string &ofname, bool transferData,
                          double sizeFactor) {
  int arrayID = IsArrayName(arrayName);
  if (arrayID == -1) {
    std::cout << "Array " << arrayName
              << " not found in set of point data arrays" << std::endl;
    exit(1);
  }
  refineMesh(method, arrayID, dev_mult, maxIsmin, edge_scale, ofname,
             transferData, sizeFactor);
}

/**
 **/
void meshBase::refineMesh(const std::string &method,
                          const std::string &arrayName, int order,
                          const std::string &ofname, bool transferData) {
  int arrayID = IsArrayName(arrayName);
  if (arrayID == -1) {
    std::cout << "Array " << arrayName
              << " not found in set of point data arrays" << std::endl;
    exit(1);
  }

  refineMesh(method, arrayID, order, ofname, transferData);
}

/** added for uniform refinement by driver
 **/
void meshBase::refineMesh(const std::string &method, double edge_scale,
                          const std::string &ofname, bool transferData,
                          bool constrainBoundary) {
  refineMesh(method, 0, 0, false, edge_scale, ofname, transferData, 1.0,
             constrainBoundary);
}

/**
 **/
vtkSmartPointer<vtkStaticCellLocator> meshBase::buildStaticCellLocator() {
  vtkSmartPointer<vtkStaticCellLocator> cellLocator =
      vtkSmartPointer<vtkStaticCellLocator>::New();
  cellLocator->SetDataSet(dataSet);
  cellLocator->BuildLocator();
  return cellLocator;
}

vtkSmartPointer<vtkStaticPointLocator> meshBase::buildStaticPointLocator() {
  auto pointLocator = vtkSmartPointer<vtkStaticPointLocator>::New();
  pointLocator->SetDataSet(dataSet);
  pointLocator->BuildLocator();
  return pointLocator;
}

/**
 **/
void meshBase::checkMesh(const std::string &ofname) const {
  std::unique_ptr<MeshQuality> qualCheck =
      std::unique_ptr<MeshQuality>(new MeshQuality(this));
  qualCheck->checkMesh(ofname);
}

/**
 **/
int diffMesh(meshBase *mesh1, meshBase *mesh2) {
  //double tol = 1e-14;
  double tol = 3e-2;

  if (mesh1->getNumberOfPoints() != mesh2->getNumberOfPoints() ||
      mesh1->getNumberOfCells() != mesh2->getNumberOfCells()) {
    std::cerr << "Meshes don't have the same number of points or cells"
              << std::endl;
    return 1;
  }

  for (int i = 0; i < mesh1->getNumberOfPoints(); ++i) {
    std::vector<double> coord1 = mesh1->getPoint(i);
    std::vector<double> coord2 = mesh2->getPoint(i);
    for (int j = 0; j < 3; ++j) {
      if (std::fabs((coord1[j] - coord2[j])/coord2[j]) > tol) {
        std::cerr << "Meshes differ in point coordinates" << std::endl;
        std::cerr << "Index " << i << " Component " << j << std::endl;
        std::cerr << "Coord 1 " << std::setprecision(15) << coord1[j]
                  << " Coord 2 " << std::setprecision(15) << coord2[j]
                  << std::endl;
        std::cerr << "Meshes differ in point coordinates" << std::endl;
        return 1;
      }
    }
  }

  for (int i = 0; i < mesh1->getNumberOfCells(); ++i) {
    std::vector<std::vector<double>> cell1 = mesh1->getCellVec(i);
    std::vector<std::vector<double>> cell2 = mesh2->getCellVec(i);
    if (cell1.size() != cell2.size()) {
      std::cerr << "Meshes differ in cells" << std::endl;
      return 1;
    }
    for (int j = 0; j < cell1.size(); ++j) {
      for (int k = 0; k < 3; ++k) {
        if (std::fabs((cell1[j][k] - cell2[j][k])/cell2[j][k]) > tol) {
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

  if (numArr1 != numArr2) {
    std::cerr << "Meshes have different numbers of point data" << std::endl;
    std::cerr << "Mesh 1 has " << numArr1 << std::endl;
    std::cerr << "Mesh 2 has " << numArr2 << std::endl;
    return 1;
  }

  for (int i = 0; i < numArr1; ++i) {
    vtkDataArray *da1 = pd1->GetArray(i);
    std::cerr << "checking array " << da1->GetName() << std::endl;
    vtkDataArray *da2 = pd2->GetArray(pd1->GetArrayName(i));
    double range[2];
    double abs_error;
    double rel_error;
    int numComponent = da1->GetNumberOfComponents();
    for (int j = 0; j < mesh1->getNumberOfPoints(); ++j) {
      double *comps1 = new double[numComponent];
      double *comps2 = new double[numComponent];
      da1->GetTuple(j, comps1);
      da2->GetTuple(j, comps2);
      for (int k = 0; k < numComponent; ++k) {
        da1->GetRange(range, k);
        abs_error = std::fabs(comps1[k] - comps2[k]);
        double max_val = std::max(std::abs(range[0]), std::abs(range[1]));
        rel_error = abs_error/std::max(1.0, max_val);
        if (rel_error > tol) {
          std::cerr << "For point data array " << da1->GetName() << std::endl;
          std::cerr << "Meshes differ in point data values at point " << j
                    << " component " << k << std::endl;
          std::cerr << std::setprecision(15)
                    << "Mesh 1 value : "
                    << comps1[k] << std::endl;
          std::cerr << std::setprecision(15)
                    << "Mesh 2 value : "
                    << comps2[k] << std::endl;
          return 1;
        }
      }
      delete[] comps1;
      delete[] comps2;
    }
  }
  std::cerr << "Meshes are the same" << std::endl;
  return 0;
}

bool sortNemId_tVec_compare::operator()(std::vector<nemId_t> lhs,
                                        std::vector<nemId_t> rhs) const {
  std::sort(lhs.begin(), lhs.end());
  std::sort(rhs.begin(), rhs.end());
  return lhs < rhs;
}

// convert all quadrilateral elements to triangular elements
// using naive quad splitting
meshBase *meshBase::convertQuads() {
  // Create new dataset
  vtkSmartPointer<vtkUnstructuredGrid> vtkDataSetNoQuads =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  // Accumulate vector of ID lists for quads removed
  std::vector<std::vector<int>> removedQuadIds;
  // Set points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(this->getDataSet()->GetNumberOfPoints());
  double pt[3];
  for (int iPoint = 0; iPoint < this->getDataSet()->GetNumberOfPoints();
       iPoint++) {
    this->getDataSet()->GetPoint(iPoint, pt);
    points->SetPoint(iPoint, pt[0], pt[1], pt[2]);
  }
  vtkDataSetNoQuads->SetPoints(points);

  // loop through elements and remove quads by splitting into tris
  for (int iCell = 0; iCell < this->getDataSet()->GetNumberOfCells(); iCell++) {
    // is quad?
    if (this->getDataSet()->GetCell(iCell)->GetCellType() == VTK_QUAD) {
      vtkSmartPointer<vtkIdList> elmIds =
          this->getDataSet()->GetCell(iCell)->GetPointIds();
      std::vector<int> elmIdsVec = {static_cast<int>(elmIds->GetId(0)),
                                    static_cast<int>(elmIds->GetId(1)),
                                    static_cast<int>(elmIds->GetId(2)),
                                    static_cast<int>(elmIds->GetId(3))};
      removedQuadIds.push_back(elmIdsVec);
    }
    // only add non-quads to vtkdataset
    else {
      vtkDataSetNoQuads->InsertNextCell(
          this->getDataSet()->GetCell(iCell)->GetCellType(),
          this->getDataSet()->GetCell(iCell)->GetPointIds());
    }
  }

  // Add in new tris based on accumulated quad Ids
  for (auto elemItr = removedQuadIds.begin(); elemItr != removedQuadIds.end();
       elemItr++) {
    vtkSmartPointer<vtkIdList> triList1 = vtkSmartPointer<vtkIdList>::New();
    triList1->InsertNextId((*elemItr)[0]);
    triList1->InsertNextId((*elemItr)[1]);
    triList1->InsertNextId((*elemItr)[2]);
    vtkDataSetNoQuads->InsertNextCell(VTK_TRIANGLE, triList1);

    vtkSmartPointer<vtkIdList> triList2 = vtkSmartPointer<vtkIdList>::New();
    triList2->InsertNextId((*elemItr)[0]);
    triList2->InsertNextId((*elemItr)[2]);
    triList2->InsertNextId((*elemItr)[3]);

    vtkDataSetNoQuads->InsertNextCell(VTK_TRIANGLE, triList2);
  }

  return Create(vtkDataSetNoQuads, "removedQuads.vtu");
}

void meshBase::convertHexToTetVTK(vtkSmartPointer<vtkDataSet> meshdataSet) {
  vtkSmartPointer<vtkDataSetTriangleFilter> triFilter =
      vtkSmartPointer<vtkDataSetTriangleFilter>::New();
  triFilter->SetInputData(meshdataSet);
  triFilter->Update();

  dataSet = triFilter->GetOutput();
}

std::vector<int> meshBase::getArrayIDs(std::vector<std::string> arrayNames,
                                       bool fromPointArrays) {
  std::vector<int> arrayIDs(arrayNames.size());
  for (int i = 0; i < arrayNames.size(); ++i) {
    int id = IsArrayName(arrayNames[i], fromPointArrays);
    if (id == -1) {
      std::cerr << "Array " << arrayNames[i]
                << " not found in set of data arrays." << std::endl;
      exit(1);
    }
    arrayIDs[i] = id;
  }
  return arrayIDs;
}
