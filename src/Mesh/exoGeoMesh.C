#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "Mesh/exoGeoMesh.H"

#include <exodusII.h>
#include <vtkCellIterator.h>
#include <vtkCompositeDataIterator.h>
#include <vtkEmptyCell.h>
#include <vtkExodusIIReader.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <numeric>
#include <unordered_set>

#include "AuxiliaryFunctions.H"
#include "MeshOperation/mergeCells.H"

#ifdef HAVE_GMSH
#include <gmsh.h>
#endif

namespace NEM {
namespace MSH {

namespace {

/**
 * Checks and prints errors, if any, from calls to the Exodus API. Note if
 * Exodus has a fatal error, will cause exit with non-zero code.
 * @param errNum return value from a call to Exodus API functions
 * @param createOrOpen whether call was from ex_create or ex_open
 */
void checkExodusErr(int errNum, bool createOrOpen = false) {
  if (errNum < 0 || (!createOrOpen && errNum > 0)) {
    const char *msg;
    const char *func;
    ex_get_err(&msg, &func, &errNum);
    ex_opts(EX_VERBOSE | EX_ABORT);
    ex_err(func, msg, errNum);
  }
}

const char *vtkCellType2exoType(VTKCellType vtkCellType) {
  switch (vtkCellType) {
    case VTK_TRIANGLE: return "TRIANGLE";
    case VTK_QUAD: return "QUAD";
    case VTK_TETRA: return "TETRA";
    case VTK_HEXAHEDRON: return "HEX";
    case VTK_WEDGE: return "WEDGE";
    default: return "OTHER";
  }
}

int exoSide2vtkSide(int exoSide, VTKCellType vtkCellType) {
  switch (vtkCellType) {
    case VTK_HEXAHEDRON:
      switch (exoSide) {
        case 1: return 2;
        case 2: return 1;
        case 3: return 3;
        case 4: return 0;
        case 5: return 4;
        case 6: return 5;
        default: return -1;
      }
    case VTK_WEDGE:
      switch (exoSide) {
        case 1: return 2;
        case 2: return 3;
        case 3: return 4;
        case 4: return 0;
        case 5: return 1;
        default: return -1;
      }
    default: return exoSide - 1;
  }
}

int vtkSide2exoSide(int vtkSide, int vtkCellType) {
  switch (vtkCellType) {
    case VTK_HEXAHEDRON:
      switch (vtkSide) {
        case 0: return 4;
        case 1: return 2;
        case 2: return 1;
        case 3: return 3;
        case 4: return 5;
        case 5: return 6;
        default: return -1;
      }
    case VTK_WEDGE:
      switch (vtkSide) {
        case 0: return 4;
        case 1: return 5;
        case 2: return 1;
        case 3: return 2;
        case 4: return 3;
        default: return -1;
      }
    default: return vtkSide + 1;
  }
}

std::string getArrComponentName(vtkDataArray *arr, int component) {
  auto componentName = arr->GetComponentName(component);
  if (componentName) {
    return std::string{componentName};
  }
  if (!arr->GetName()) {
    return std::string{};
  }
  std::string varName = arr->GetName();
  if (arr->GetNumberOfComponents() == 2) {
    if (component == 0) {
      varName += 'r';
    } else {  // component == 1
      varName += 'z';
    }
  } else if (arr->GetNumberOfComponents() == 3) {
    if (component == 0) {
      varName += 'x';
    } else if (component == 1) {
      varName += 'y';
    } else {  // component == 2
      varName += 'z';
    }
  } else if (arr->GetNumberOfComponents() == 6) {
    if (component == 0) {
      varName += "xx";
    } else if (component == 1) {
      varName += "yy";
    } else if (component == 2) {
      varName += "zz";
    } else if (component == 3) {
      varName += "xy";
    } else if (component == 4) {
      varName += "xz";
    } else {  // component == 5
      varName += "yz";
    }
  } else if (arr->GetNumberOfComponents() > 1) {
    varName += std::to_string(component);
  }
  return varName;
}

/**
 * Combine @p mesh1 and @p mesh2 into @p outMesh. PointData and CellData are
 * added to @p outMesh only if present in both @p mesh1 and @p mesh2.
 * @param mesh1 First mesh to combine. Note that cell global IDs and cell
 * pedigree IDs are removed, if present.
 * @param mesh2 Second mesh to combine.
 * @param outMesh Resulting vtkUnstructuredGrid. Should have no points, cells,
 * vtkPointData, or vtkCellData.
 * @param tol tolerance for merging points
 * @return maps from mesh1 points to outMesh points and mesh2 points to outMesh
 * points (all using 0-based indices; point global IDs are not used)
 */
std::pair<vtkSmartPointer<vtkIdTypeArray>, vtkSmartPointer<vtkIdTypeArray>>
mergeMeshes(vtkUnstructuredGrid *mesh1, vtkUnstructuredGrid *mesh2,
            vtkUnstructuredGrid *outMesh, float tol) {
  if (mesh1->GetCellData()->GetGlobalIds()) {
    mesh1->GetCellData()->RemoveArray(
        mesh1->GetCellData()->GetGlobalIds()->GetName());
  }
  if (mesh1->GetCellData()->GetPedigreeIds()) {
    mesh1->GetCellData()->RemoveArray(
        mesh1->GetCellData()->GetPedigreeIds()->GetName());
  }
  // Note that mesh1's Points will be overwritten in first call to MergeDataSet
  auto merge = vtkSmartPointer<mergeCells>::New();
  merge->SetUnstructuredGrid(outMesh);
  merge->SetTotalNumberOfDataSets(2);
  merge->SetTotalNumberOfPoints(mesh1->GetNumberOfPoints() +
                                mesh2->GetNumberOfPoints());
  merge->SetTotalNumberOfCells(mesh1->GetNumberOfCells() +
                               mesh2->GetNumberOfCells());
  merge->SetPointMergeTolerance(tol);
  auto nodeMap1 = merge->MergeDataSet(mesh1);
  auto nodeMap2 = merge->MergeDataSet(mesh2);
  merge->Finish();
  for (int i = outMesh->GetCellData()->GetNumberOfArrays() - 1; i >= 0; --i) {
    if (outMesh->GetCellData()->GetArray(i)->GetNumberOfTuples() !=
        outMesh->GetNumberOfCells()) {
      outMesh->GetCellData()->RemoveArray(i);
    }
  }
  for (int i = outMesh->GetPointData()->GetNumberOfArrays() - 1; i >= 0; --i) {
    if (outMesh->GetPointData()->GetArray(i)->GetNumberOfTuples() !=
        outMesh->GetNumberOfPoints()) {
      outMesh->GetPointData()->RemoveArray(i);
    }
  }
  return {nodeMap1, nodeMap2};
}

void resetSideSetPoints(vtkUnstructuredGrid *mesh, vtkPolyData *sideSet,
                        vtkIdTypeArray *nodeMap) {
  if (sideSet) {
    sideSet->SetPoints(mesh->GetPoints());
    for (vtkIdType i = 0; i < sideSet->GetNumberOfCells(); ++i) {
      auto oldCell = sideSet->GetCell(i);
      auto numPoints = oldCell->GetNumberOfPoints();
      auto oldPoints = oldCell->GetPointIds();
      for (int j = 0; j < numPoints; ++j) {
        oldPoints->SetId(j, nodeMap->GetValue(oldPoints->GetId(j)));
      }
      sideSet->ReplaceCell(i, numPoints, oldPoints->GetPointer(0));
    }
  }
}

}  // namespace

vtkStandardNewMacro(exoGeoMesh)

exoGeoMesh *exoGeoMesh::Read(const std::string &fileName, int timeStep) {
  return new exoGeoMesh(fileName, timeStep);
}

const char *exoGeoMesh::getExoIdArrName() {
  return vtkExodusIIReader::GetObjectIdArrayName();
}

exoGeoMesh::exoGeoMesh() : exoGeoMesh("") {}

exoGeoMesh::exoGeoMesh(const vtkSmartPointer<vtkExodusIIReader> &reader)
    : geoMeshBase(exoReader2GM(reader)) {
  if (reader && reader->GetOutput() &&
      reader->GetOutput()->GetNumberOfBlocks() == 8) {
    auto mbds = reader->GetOutput();
    _title = reader->GetTitle();
    auto elemBlocks = vtkMultiBlockDataSet::SafeDownCast(mbds->GetBlock(0));
    for (unsigned int i = 0; i < elemBlocks->GetNumberOfBlocks(); ++i) {
      std::string elemBlockName =
          elemBlocks->GetMetaData(i)->Get(vtkMultiBlockDataSet::NAME());
      auto elemBlock =
          vtkUnstructuredGrid::SafeDownCast(elemBlocks->GetBlock(i));
      int elemBlockId;
      VTKCellType cellType;
      if (elemBlock->GetNumberOfCells() > 0) {
        elemBlockId =
            vtkIntArray::FastDownCast(
                elemBlock->GetCellData()->GetAbstractArray(getExoIdArrName()))
                ->GetValue(0);
        cellType = static_cast<VTKCellType>(elemBlock->GetCellType(0));
      } else {
        elemBlockId =
            vtkIntArray::FastDownCast(
                elemBlock->GetFieldData()->GetAbstractArray("ElementBlockIds"))
                ->GetValue(0);
        cellType = VTK_EMPTY_CELL;
      }
      _elemBlocks[elemBlockId] = {elemBlockName, cellType, {}};
    }
    auto sideSets = vtkMultiBlockDataSet::SafeDownCast(mbds->GetBlock(4));
    if (sideSets->GetNumberOfBlocks() >= 1) {
      for (int i = 0; i < sideSets->GetNumberOfBlocks(); ++i) {
        auto sideSetName = reader->GetSideSetArrayName(i);
        auto sideSetId = reader->GetObjectId(vtkExodusIIReader::SIDE_SET, i);
        _sideSetNames[sideSetId] = sideSetName ? sideSetName : "";
      }
    }

    auto mesh = getGeoMesh().mesh;
    std::map<vtkIdType, vtkIdType> globalNodeId2meshIdx;
    auto globalNodeIdArr =
        vtkIdTypeArray::FastDownCast(mesh->GetPointData()->GetAbstractArray(
            reader->GetGlobalNodeIdArrayName()));
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      globalNodeId2meshIdx[globalNodeIdArr->GetValue(i)] = i;
    }
    mesh->GetPointData()->RemoveArray(reader->GetGlobalNodeIdArrayName());
    mesh->GetPointData()->RemoveArray(reader->GetPedigreeNodeIdArrayName());
    auto nodeSets = vtkMultiBlockDataSet::SafeDownCast(mbds->GetBlock(7));
    for (unsigned int i = 0; i < nodeSets->GetNumberOfBlocks(); ++i) {
      std::string nodeSetName =
          nodeSets->GetMetaData(i)->Get(vtkMultiBlockDataSet::NAME());
      auto nodeSet = vtkUnstructuredGrid::SafeDownCast(nodeSets->GetBlock(i));
      auto nodeSetId =
          vtkIntArray::FastDownCast(
              nodeSet->GetCellData()->GetAbstractArray(getExoIdArrName()))
              ->GetValue(0);
      std::vector<vtkIdType> nodes;
      nodes.reserve(nodeSet->GetNumberOfCells());
      auto nodeSetGlobalNodeIdArr = vtkIdTypeArray::FastDownCast(
          nodeSet->GetPointData()->GetAbstractArray(
              reader->GetGlobalNodeIdArrayName()));
      for (vtkIdType j = 0; j < nodeSet->GetNumberOfCells(); ++j) {
        auto nodeGlobalId = nodeSetGlobalNodeIdArr->GetValue(
            nodeSet->GetCell(j)->GetPointId(0));
        auto find = globalNodeId2meshIdx.find(nodeGlobalId);
        if (find != globalNodeId2meshIdx.end()) {
          nodes.emplace_back(find->second);
        }
      }
      _nodeSets[nodeSetId] = {nodeSetName, nodes};
    }
    // Unfortunately the vtkExodusIIReader does not seem to support exodus
    // properties at the moment.
    int comp_ws = sizeof(double);
    int io_ws = 0;
    float version;
    auto exoid =
        ex_open(reader->GetFileName(), EX_READ, &comp_ws, &io_ws, &version);
    checkExodusErr(exoid, true);
    // We get the properties in order, so, just in case, we get the element
    // block ids in order
    std::vector<int> elemBlockIds;
    elemBlockIds.resize(_elemBlocks.size());
    auto err = ex_get_ids(exoid, EX_ELEM_BLOCK, elemBlockIds.data());
    checkExodusErr(err);
    int numProps = ex_inquire_int(exoid, EX_INQ_EB_PROP);
    std::vector<std::string> propNames;
    propNames.resize(numProps);
    std::vector<const char *> propNamesPtr;
    propNamesPtr.reserve(numProps);
    for (int i = 0; i < numProps; ++i) {
      propNames[i].resize(MAX_NAME_LENGTH);
      propNamesPtr.emplace_back(propNames[i].c_str());
    }
    err = ex_get_prop_names(exoid, EX_ELEM_BLOCK,
                            const_cast<char **>(propNamesPtr.data()));
    checkExodusErr(err);
    for (int i = 0; i < numProps; ++i) {
      if (strcmp(propNamesPtr[i], "ID") != 0) {
        std::vector<int> propVals;
        propVals.resize(_elemBlocks.size());
        err = ex_get_prop_array(exoid, EX_ELEM_BLOCK, propNamesPtr[i],
                                propVals.data());
        checkExodusErr(err);
        _elemBlockPropNames.insert(propNamesPtr[i]);
        for (std::size_t j = 0; j < _elemBlocks.size(); ++j) {
          _elemBlocks[elemBlockIds[j]].properties[propNamesPtr[i]] =
              propVals[j];
        }
      }
    }
    err = ex_close(exoid);
    checkExodusErr(err);
  }
  std::cout << "exoGeoMesh constructed" << std::endl;
}

exoGeoMesh::exoGeoMesh(const std::string &fileName, int timeStep)
    : exoGeoMesh(getExoReader(fileName, timeStep)) {}

exoGeoMesh::~exoGeoMesh() { std::cout << "exoGeoMesh destructed" << std::endl; }

void exoGeoMesh::write(const std::string &fileName) {
  std::string fileExt = nemAux::find_ext(fileName);
  if (fileExt == ".exo" || fileExt == ".e" || fileExt == ".gen" ||
      fileExt == ".g") {
    int comp_ws = sizeof(double);
    int io_ws = sizeof(double);
    int exoid = ex_create(fileName.c_str(), EX_CLOBBER, &comp_ws, &io_ws);
    checkExodusErr(exoid, true);

    auto geoMesh = this->getGeoMesh();
    auto mesh = geoMesh.mesh;
    auto sideSet = geoMesh.sideSet;
    if (!mesh || mesh->GetNumberOfCells() == 0) {
      std::cerr << "No mesh information to write." << std::endl;
      return;
    }

    // Make sure no empty element blocks
    std::map<int, std::vector<vtkIdType>> elemBlocks;
    {
      auto meshEntities = vtkIntArray::FastDownCast(
          mesh->GetCellData()->GetAbstractArray(getExoIdArrName()));
      for (const auto &elemBlock : _elemBlocks) {
        elemBlocks[elemBlock.first];
      }
      for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
        auto entity = meshEntities->GetValue(i);
        elemBlocks.at(entity).emplace_back(i);
      }
      for (auto it = elemBlocks.begin(); it != elemBlocks.end();) {
        if (it->second.empty()) {
          it = elemBlocks.erase(it);
        } else {
          ++it;
        }
      }
    }

    // Get number of element block variables, nodal variables, and global
    // variables
    int numElemVars = 0;
    int numNodeVars = 0;
    int numGlobalVars = 0;
    std::vector<bool> cellDataWriteFlag;
    std::vector<bool> pointDataWriteFlag;
    std::vector<bool> fieldDataWriteFlag;
    int timeStep = 1;
    {
      // Cell Data
      auto globalIds = mesh->GetCellData()->GetGlobalIds();
      const char *globalIdName;
      if (globalIds) {
        globalIdName = globalIds->GetName();
      } else {
        globalIdName = "GlobalElementId";
      }
      auto pedigreeIds = mesh->GetCellData()->GetPedigreeIds();
      const char *pedgreeIdName;
      if (pedigreeIds) {
        pedgreeIdName = pedigreeIds->GetName();
      } else {
        pedgreeIdName = "PedigreeElementId";
      }
      auto numCells = mesh->GetNumberOfCells();
      cellDataWriteFlag.resize(mesh->GetCellData()->GetNumberOfArrays(), false);
      for (int i = 0; i < mesh->GetCellData()->GetNumberOfArrays(); ++i) {
        // Use GetArray instead of GetAbstractArray to filter for vtkDataArray
        auto arr = mesh->GetCellData()->GetArray(i);
        if (arr) {
          auto arrName = arr->GetName();
          if (strcmp(arrName, globalIdName) != 0 &&
              strcmp(arrName, pedgreeIdName) != 0 &&
              strcmp(arrName, getExoIdArrName()) != 0 &&
              arr->GetNumberOfTuples() == numCells) {
            auto numComponents = arr->GetNumberOfComponents();
            cellDataWriteFlag[i] = (numComponents > 0);
            numElemVars += numComponents;
          }
        }
      }
      // Point Data
      globalIds = mesh->GetPointData()->GetGlobalIds();
      if (globalIds) {
        globalIdName = globalIds->GetName();
      } else {
        globalIdName = "GlobalNodeId";
      }
      pedigreeIds = mesh->GetPointData()->GetPedigreeIds();
      if (pedigreeIds) {
        pedgreeIdName = pedigreeIds->GetName();
      } else {
        pedgreeIdName = "PedigreeNodeId";
      }
      auto numNodes = mesh->GetNumberOfPoints();
      pointDataWriteFlag.resize(mesh->GetPointData()->GetNumberOfArrays(),
                                false);
      for (int i = 0; i < mesh->GetPointData()->GetNumberOfArrays(); ++i) {
        // Use GetArray instead of GetAbstractArray to filter for vtkDataArray
        auto arr = mesh->GetPointData()->GetArray(i);
        if (arr) {
          auto arrName = arr->GetName();
          if (strcmp(arrName, globalIdName) != 0 &&
              strcmp(arrName, pedgreeIdName) != 0 &&
              arr->GetNumberOfTuples() == numNodes) {
            auto numComponents = arr->GetNumberOfComponents();
            pointDataWriteFlag[i] = (numComponents > 0);
            numNodeVars += numComponents;
          }
        }
      }
      // Field Data
      fieldDataWriteFlag.resize(mesh->GetFieldData()->GetNumberOfArrays(),
                                false);
      for (int i = 0; i < mesh->GetFieldData()->GetNumberOfArrays(); ++i) {
        // Use GetArray instead of GetAbstractArray to filter for vtkDataArray
        auto arr = mesh->GetFieldData()->GetArray(i);
        if (arr && arr->GetNumberOfTuples() == 1) {
          auto numComponents = arr->GetNumberOfComponents();
          fieldDataWriteFlag[i] = (numComponents > 0);
          numGlobalVars += numComponents;
        }
      }
    }

    // Have to make sure no empty side sets; might as well get the
    // indices now.
    std::map<int, std::vector<vtkIdType>> exoSideSets;
    if (sideSet.sides) {
      auto sideSetExoId = vtkIntArray::FastDownCast(
          sideSet.sides->GetCellData()->GetAbstractArray(getExoIdArrName()));
      for (vtkIdType i = 0; i < sideSet.sides->GetNumberOfCells(); ++i) {
        auto entity = sideSetExoId->GetValue(i);
        exoSideSets[entity].emplace_back(i);
      }
    }

    // Make sure no empty node sets and no duplicate nodes
    resetNodeSetPoints();
    int numNodeSets = std::accumulate(
        _nodeSets.begin(), _nodeSets.end(), 0,
        [](const int &a, const std::pair<const int, nodeSet> &b) {
          return a + (!b.second.nodes.empty());
        });

    auto err = ex_put_init(exoid, _title.c_str(), geoMesh.getDimension(),
                           mesh->GetNumberOfPoints(), mesh->GetNumberOfCells(),
                           elemBlocks.size(), numNodeSets, exoSideSets.size());
    checkExodusErr(err);
    if (numElemVars > 0) {
      err = ex_put_variable_param(exoid, EX_ELEM_BLOCK, numElemVars);
      checkExodusErr(err);
    }
    if (numNodeVars > 0) {
      err = ex_put_variable_param(exoid, EX_NODAL, numNodeVars);
      checkExodusErr(err);
    }
    if (numGlobalVars > 0) {
      err = ex_put_variable_param(exoid, EX_GLOBAL, numGlobalVars);
      checkExodusErr(err);
    }
    if (numElemVars > 0 || numNodeVars > 0 || numGlobalVars > 0) {
      double time = 0;
      err = ex_put_time(exoid, timeStep, &time);
      checkExodusErr(err);
    }

    std::vector<double> x_coord;
    std::vector<double> y_coord;
    std::vector<double> z_coord;
    x_coord.resize(mesh->GetNumberOfPoints(), 0);
    y_coord.resize(mesh->GetNumberOfPoints(), 0);
    z_coord.resize(mesh->GetNumberOfPoints(), 0);
    auto coord = mesh->GetPoints()->GetData();
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      x_coord[i] = coord->GetComponent(i, 0);
      y_coord[i] = coord->GetComponent(i, 1);
      z_coord[i] = coord->GetComponent(i, 2);
    }
    err = ex_put_coord(exoid, x_coord.data(), y_coord.data(), z_coord.data());
    checkExodusErr(err);

    // map from index in mesh to implicit index in exodus output
    std::vector<vtkIdType> vtkIdMap(mesh->GetNumberOfCells(), -1);
    int lastElemIdx = 0;
    for (const auto &elemBlock : elemBlocks) {
      auto firstCell = mesh->GetCell(elemBlock.second[0]);
      auto cellType = static_cast<VTKCellType>(firstCell->GetCellType());
      err = ex_put_block(exoid, EX_ELEM_BLOCK, elemBlock.first,
                         vtkCellType2exoType(cellType), elemBlock.second.size(),
                         firstCell->GetNumberOfPoints(), 0, 0, 0);
      checkExodusErr(err);
      std::vector<int> conn;
      conn.reserve(elemBlock.second.size() * firstCell->GetNumberOfPoints());
      for (auto elem : elemBlock.second) {
        auto cell = mesh->GetCell(elem);
        vtkIdMap[elem] = ++lastElemIdx;
        for (int i = 0; i < cell->GetNumberOfPoints(); ++i) {
          // Exodus starts indexing from 1
          conn.emplace_back(cell->GetPointId(i) + 1);
        }
      }
      err = ex_put_conn(exoid, EX_ELEM_BLOCK, elemBlock.first, conn.data(),
                        nullptr, nullptr);
      checkExodusErr(err);
      if (!_elemBlocks.at(elemBlock.first).name.empty()) {
        err = ex_put_name(exoid, EX_ELEM_BLOCK, elemBlock.first,
                          _elemBlocks.at(elemBlock.first).name.c_str());
        checkExodusErr(err);
      }
    }

    if (sideSet.sides) {
      auto sideSetOrigCellId = geoMesh.sideSet.getOrigCellArr();
      auto sideSetCellFaceId = geoMesh.sideSet.getCellFaceArr();
      for (const auto &exoSideSet : exoSideSets) {
        err = ex_put_set_param(exoid, EX_SIDE_SET, exoSideSet.first,
                               exoSideSet.second.size(), 0);
        checkExodusErr(err);
        std::vector<int> sideSetElem, sideSetSide;
        sideSetElem.reserve(exoSideSet.second.size());
        sideSetSide.reserve(exoSideSet.second.size());
        for (auto sideIdx : exoSideSet.second) {
          auto parentCellIdx = sideSetOrigCellId->GetTypedComponent(sideIdx, 0);
          sideSetElem.emplace_back(vtkIdMap[parentCellIdx]);
          sideSetSide.emplace_back(
              vtkSide2exoSide(sideSetCellFaceId->GetTypedComponent(sideIdx, 0),
                              mesh->GetCell(parentCellIdx)->GetCellType()));
        }
        err = ex_put_set(exoid, EX_SIDE_SET, exoSideSet.first,
                         sideSetElem.data(), sideSetSide.data());
        checkExodusErr(err);
        auto foundName = _sideSetNames.find(exoSideSet.first);
        if (foundName != _sideSetNames.end()) {
          err = ex_put_name(exoid, EX_SIDE_SET, foundName->first,
                            foundName->second.c_str());
          checkExodusErr(err);
        }
      }
    }

    for (const auto &nodeSet : _nodeSets) {
      if (!nodeSet.second.nodes.empty()) {
        err = ex_put_set_param(exoid, EX_NODE_SET, nodeSet.first,
                               nodeSet.second.nodes.size(), 0);
        checkExodusErr(err);
        std::vector<int> nodes;  // Have to add 1
        nodes.reserve(nodeSet.second.nodes.size());
        std::transform(nodeSet.second.nodes.begin(), nodeSet.second.nodes.end(),
                       std::back_inserter(nodes), [](int x) { return x + 1; });
        err = ex_put_set(exoid, EX_NODE_SET, nodeSet.first, nodes.data(),
                         nullptr);
        checkExodusErr(err);
        if (!nodeSet.second.name.empty()) {
          err = ex_put_name(exoid, EX_NODE_SET, nodeSet.first,
                            nodeSet.second.name.c_str());
          checkExodusErr(err);
        }
      }
    }

    if (numElemVars > 0) {
      int var_index = 1;
      std::vector<double> vals(mesh->GetNumberOfCells());
      for (int i = 0; i < mesh->GetCellData()->GetNumberOfArrays(); ++i) {
        if (cellDataWriteFlag[i]) {
          auto arr = mesh->GetCellData()->GetArray(i);
          for (int j = 0; j < arr->GetNumberOfComponents(); ++j) {
            for (vtkIdType k = 0; k < mesh->GetNumberOfCells(); ++k) {
              vals[vtkIdMap[k] - 1] = arr->GetComponent(k, j);
            }
            auto name = getArrComponentName(arr, j);
            if (name.empty()) {
              name = "ElemBlockVar" + std::to_string(var_index);
            }
            err = ex_put_variable_name(exoid, EX_ELEM_BLOCK, var_index,
                                       name.c_str());
            checkExodusErr(err);
            std::size_t varOffset = 0;
            for (const auto &elemBlock : elemBlocks) {
              err = ex_put_var(exoid, timeStep, EX_ELEM_BLOCK, var_index,
                               elemBlock.first, elemBlock.second.size(),
                               vals.data() + varOffset);
              checkExodusErr(err);
              varOffset += elemBlock.second.size();
            }
            ++var_index;
          }
        }
      }
    }

    if (numNodeVars > 0) {
      int var_index = 1;
      std::vector<double> vals(mesh->GetNumberOfPoints());
      for (int i = 0; i < mesh->GetPointData()->GetNumberOfArrays(); ++i) {
        if (pointDataWriteFlag[i]) {
          auto arr = mesh->GetPointData()->GetArray(i);
          for (int j = 0; j < arr->GetNumberOfComponents(); ++j) {
            auto name = getArrComponentName(arr, j);
            if (name.empty()) {
              name = "NodalVar" + std::to_string(var_index);
            }
            err =
                ex_put_variable_name(exoid, EX_NODAL, var_index, name.c_str());
            checkExodusErr(err);
            for (vtkIdType k = 0; k < mesh->GetNumberOfPoints(); ++k) {
              vals[k] = arr->GetComponent(k, j);
            }
            err = ex_put_var(exoid, timeStep, EX_NODAL, var_index, 1,
                             mesh->GetNumberOfPoints(), vals.data());
            checkExodusErr(err);
            ++var_index;
          }
        }
      }
    }
    if (numGlobalVars > 0) {
      int var_index = 1;
      std::vector<double> vals;
      for (int i = 0; i < mesh->GetFieldData()->GetNumberOfArrays(); ++i) {
        if (fieldDataWriteFlag[i]) {
          auto arr = mesh->GetFieldData()->GetArray(i);
          for (int j = 0; j < arr->GetNumberOfComponents(); ++j) {
            vals.emplace_back(arr->GetComponent(0, j));
            auto name = getArrComponentName(arr, j);
            if (!name.empty()) {
              err = ex_put_variable_name(exoid, EX_GLOBAL, var_index,
                                         name.c_str());
              checkExodusErr(err);
            }
            ++var_index;
          }
        }
      }
      err = ex_put_var(exoid, timeStep, EX_GLOBAL, 1, 0, numGlobalVars,
                       vals.data());
      checkExodusErr(err);
    }

    if (!_elemBlockPropNames.empty()) {
      // Element block properties
      std::vector<const char *> propNames;
      propNames.reserve(_elemBlockPropNames.size());
      std::vector<std::vector<int>> propVals;
      propVals.resize(_elemBlockPropNames.size());
      for (const auto &elemBlockProp : _elemBlockPropNames) {
        propNames.emplace_back(elemBlockProp.c_str());
      }
      err = ex_put_prop_names(exoid, EX_ELEM_BLOCK, _elemBlockPropNames.size(),
                              const_cast<char **>(propNames.data()));
      checkExodusErr(err);
      // Use elemBlocks to verify non-empty
      auto nonEmptyIt = elemBlocks.begin();
      auto it = _elemBlocks.begin();
      while (nonEmptyIt != elemBlocks.end()) {
        if (it->first == nonEmptyIt->first) {
          std::size_t propValsIdx = 0;
          auto propNamesIt = propNames.begin();
          auto propValsIt = it->second.properties.begin();
          while (propNamesIt != propNames.end()) {
            if (propValsIt != it->second.properties.end() &&
                *propNamesIt == propValsIt->first) {
              propVals[propValsIdx].emplace_back(propValsIt->second);
              ++propValsIt;
            } else {
              propVals[propValsIdx].emplace_back(0);
            }
            ++propNamesIt;
            ++propValsIdx;
          }
          ++nonEmptyIt;
        }
        ++it;
      }
      for (std::size_t i = 0; i < propNames.size(); ++i) {
        err = ex_put_prop_array(exoid, EX_ELEM_BLOCK, propNames[i],
                                propVals[i].data());
        checkExodusErr(err);
      }
    }
    err = ex_close(exoid);
    checkExodusErr(err);
  } else {
    std::cerr << "For exoGeoMesh, " << fileExt
              << " format not currently supported." << std::endl;
  }
}

void exoGeoMesh::report(std::ostream &out) const {
  geoMeshBase::report(out);
  out << "Title:\t" << _title << '\n';
  out << "Element Blocks:\t" << _elemBlocks.size() << '\n';
  out << "Side Sets:\t" << _sideSetNames.size() << '\n';
  out << "Node Sets:\t" << _nodeSets.size() << '\n';
}

void exoGeoMesh::takeGeoMesh(geoMeshBase *otherGeoMesh) {
  auto otherExoGM = dynamic_cast<exoGeoMesh *>(otherGeoMesh);
  if (otherExoGM) {
    setGeoMesh(otherExoGM->getGeoMesh());
    otherExoGM->setGeoMesh(
        {vtkSmartPointer<vtkUnstructuredGrid>::New(), {}, {}, {}});
    _title = std::move(otherExoGM->_title);
    _elemBlocks = std::move(otherExoGM->_elemBlocks);
    _sideSetNames = std::move(otherExoGM->_sideSetNames);
    _nodeSets = std::move(otherExoGM->_nodeSets);
    _elemBlockPropNames = std::move(otherExoGM->_elemBlockPropNames);
    _physGrpName = std::move(otherExoGM->_physGrpName);
    otherExoGM->resetNative();
  } else {
    geoMeshBase::takeGeoMesh(otherGeoMesh);
  }
}

void exoGeoMesh::reconstructGeo() {
  std::string link =
      getGeoMesh().link.empty() ? GEO_ENT_DEFAULT_NAME : getGeoMesh().link;
  auto mesh = getGeoMesh().mesh;
  vtkSmartPointer<vtkDataArray> linkArr =
      mesh->GetCellData()->GetArray(link.c_str());
  if (!linkArr) {
    linkArr = vtkSmartPointer<vtkIntArray>::New();
    linkArr->SetName(link.c_str());
    linkArr->SetNumberOfTuples(mesh->GetNumberOfCells());
    mesh->GetCellData()->AddArray(linkArr);
  }
  auto elemBlockIdArr = vtkIntArray::FastDownCast(
      mesh->GetCellData()->GetAbstractArray(getExoIdArrName()));
  for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
    auto physGroup = elemBlockIdArr->GetValue(i);
    if (!_physGrpName.empty()) {
      physGroup = getElemBlockProperty(_physGrpName, physGroup);
    }
    linkArr->SetComponent(i, 0, physGroup);
  }
  setGeoMesh({getGeoMesh().mesh, getGeoMesh().geo, link, getGeoMesh().sideSet});
  geoMeshBase::reconstructGeo();
  // Should now always have side set
  setSideSetObjId();
}

const std::string &exoGeoMesh::getTitle() const { return _title; }

void exoGeoMesh::setTitle(const std::string &title) { _title = title; }

int exoGeoMesh::getNumElemBlocks() const { return _elemBlocks.size(); }

std::map<int, std::string> exoGeoMesh::getElemBlockNames() const {
  std::map<int, std::string> names;
  std::transform(_elemBlocks.begin(), _elemBlocks.end(),
                 std::inserter(names, names.begin()),
                 [](const decltype(_elemBlocks)::value_type &block)
                     -> decltype(names)::value_type {
                   return {block.first, block.second.name};
                 });
  return names;
}

const std::string &exoGeoMesh::getElemBlockName(int id) const {
  auto it = _elemBlocks.find(id);
  if (it == _elemBlocks.end()) {
    std::cerr << "No block found with id: " << id << std::endl;
    exit(1);
  } else {
    return it->second.name;
  }
}

void exoGeoMesh::setElemBlockName(int id, const std::string &name) {
  auto it = _elemBlocks.find(id);
  if (it == _elemBlocks.end()) {
    std::cerr << "No block found with id: " << id << std::endl;
    exit(1);
  } else {
    it->second.name = name;
  }
}

std::vector<int> exoGeoMesh::getElemBlockIds() const {
  std::vector<int> ids;
  ids.reserve(_elemBlocks.size());
  for (const auto &elemBlock : _elemBlocks) {
    ids.emplace_back(elemBlock.first);
  }
  return ids;
}

std::vector<vtkIdType> exoGeoMesh::getElemBlock(int blockId) const {
  auto blockIt = _elemBlocks.find(blockId);
  if (blockIt == _elemBlocks.end()) {
    std::cerr << "No block found with id: " << blockId << std::endl;
    exit(1);
  }
  std::vector<vtkIdType> cellIdxList;
  auto elemBlockIdArr = vtkIntArray::FastDownCast(
      getGeoMesh().mesh->GetCellData()->GetAbstractArray(getExoIdArrName()));
  for (vtkIdType i = 0; i < getGeoMesh().mesh->GetNumberOfCells(); ++i) {
    if (elemBlockIdArr->GetValue(i) == blockId) {
      cellIdxList.emplace_back(i);
    }
  }
  return cellIdxList;
}

int exoGeoMesh::getElemBlockId(vtkIdType cellIdx) const {
  auto mesh = getGeoMesh().mesh;
  if (!mesh || mesh->GetNumberOfCells() == 0) {
    std::cerr << "Mesh is empty. No cells to find" << std::endl;
    return -1;
  }
  if (cellIdx < 0 || cellIdx >= mesh->GetNumberOfCells()) {
    std::cerr << "No cell found with index: " << cellIdx << std::endl;
    return -1;
  }
  auto elemBlockIdArr = vtkIntArray::FastDownCast(
      mesh->GetCellData()->GetAbstractArray(getExoIdArrName()));
  return elemBlockIdArr->GetValue(cellIdx);
}

int exoGeoMesh::addElemBlock(vtkUnstructuredGrid *elements,
                             const std::string &name, float tol) {
  if (elements && elements->GetNumberOfCells() > 0) {
    auto types = vtkSmartPointer<vtkCellTypes>::New();
    elements->GetCellTypes(types);
    if (types->GetNumberOfTypes() > 1) {
      std::cerr << "Single element block cannot have multiple cell types"
                << std::endl;
      return -1;
    }
    auto mesh = getGeoMesh().mesh;
    if (mesh->GetNumberOfCells() == 0) {
      mesh->DeepCopy(elements);
#ifdef HAVE_GMSH
      if (!getGeoMesh().geo.empty()) {
        gmsh::model::setCurrent(getGeoMesh().geo);
        gmsh::model::remove();
      }
#endif
      setGeoMesh({mesh, "", "", {}});
      resetNative();
      // Assumption that there is exactly one element block now.
      return _elemBlocks.begin()->first;
    } else {
      int newElemBlockId = nemAux::leastUnusedKey(_elemBlocks);
      auto &elemBlock = _elemBlocks[newElemBlockId];
      elemBlock.cellType = static_cast<VTKCellType>(types->GetCellType(0));
      elemBlock.name = name;
      addCellsToBlock(elements, newElemBlockId, tol);
      return newElemBlockId;
    }
  } else {
    std::cerr << "No elements to add." << std::endl;
    return -1;
  }
}

void exoGeoMesh::addCellsToBlock(vtkUnstructuredGrid *elements, int blockId,
                                 float tol) {
  if (elements && elements->GetNumberOfCells() > 0) {
    auto types = vtkSmartPointer<vtkCellTypes>::New();
    elements->GetCellTypes(types);
    if (types->GetNumberOfTypes() > 1) {
      std::cerr << "Single element block cannot have multiple cell types"
                << std::endl;
      exit(1);
    }
    auto elemBlock = _elemBlocks.find(blockId);
    if (getGeoMesh().mesh->GetNumberOfCells() == 0) {
      _elemBlocks.clear();
      std::cerr << "No existing element blocks." << std::endl;
      exit(1);
    }
    if (elemBlock == _elemBlocks.end()) {
      std::cerr << "No existing element block with id: " << blockId
                << std::endl;
      exit(1);
    }
    if (types->GetCellType(0) != elemBlock->second.cellType) {
      std::cerr
          << "Cell type does not match cell type of existing element block."
          << std::endl;
      exit(1);
    }
    vtkSmartPointer<vtkAbstractArray> oldBlockIdArr =
        elements->GetCellData()->GetAbstractArray(getExoIdArrName());
    if (oldBlockIdArr) {
      elements->GetCellData()->RemoveArray(getExoIdArrName());
    }
    auto mesh = getGeoMesh().mesh;
    auto newBlockIdArr = vtkSmartPointer<vtkIntArray>::Take(
        vtkIntArray::FastDownCast(mesh->GetCellData()
                                      ->GetAbstractArray(getExoIdArrName())
                                      ->NewInstance()));
    newBlockIdArr->SetName(getExoIdArrName());
    newBlockIdArr->SetNumberOfValues(elements->GetNumberOfCells());
    newBlockIdArr->FillComponent(0, blockId);
    elements->GetCellData()->AddArray(newBlockIdArr);
#ifdef HAVE_GMSH
    if (!getGeoMesh().geo.empty()) {
      gmsh::model::setCurrent(getGeoMesh().geo);
      gmsh::model::remove();
      setGeoMesh({mesh, "", "", getGeoMesh().sideSet});
    }
#endif
    auto newMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto nodeMaps = mergeMeshes(mesh, elements, newMesh, tol);
    setGeoMesh(
        {newMesh, getGeoMesh().geo, getGeoMesh().link, getGeoMesh().sideSet});
    resetNodeSetPoints(nodeMaps.first);
    resetSideSetPoints(newMesh, getGeoMesh().sideSet.sides, nodeMaps.first);
    elements->GetCellData()->RemoveArray(getExoIdArrName());
    if (oldBlockIdArr) {
      elements->GetCellData()->AddArray(oldBlockIdArr);
    }
  } else {
    std::cerr << "No elements to add." << std::endl;
  }
}

int exoGeoMesh::reassignCells(const std::vector<vtkIdType> &cells,
                              int blockId, const std::string &elemBlockName) {
  if (!cells.empty()) {
    auto mesh = getGeoMesh().mesh;
    auto cellType =
        static_cast<VTKCellType>(mesh->GetCellType(cells.front()));
    for (const auto &cell : cells) {
      if (mesh->GetCellType(cell) != cellType) {
        std::cerr << "Single element block cannot have multiple cell types"
                  << std::endl;
        exit(1);
      }
    }
    auto elemBlockIt = _elemBlocks.find(blockId);
    if (elemBlockIt == _elemBlocks.end()) {
      if (blockId <= 0) {
        blockId = nemAux::leastUnusedKey(_elemBlocks);
      }
      _elemBlocks[blockId] = {elemBlockName, cellType, {}};
    }
    auto elemBlockIdArr = vtkIntArray::FastDownCast(
        mesh->GetCellData()->GetAbstractArray(getExoIdArrName()));
    for (const auto &cell : cells) {
      elemBlockIdArr->SetValue(cell, blockId);
    }
#ifdef HAVE_GMSH
    if (!getGeoMesh().geo.empty()) {
      gmsh::model::setCurrent(getGeoMesh().geo);
      gmsh::model::remove();
      setGeoMesh({mesh, "", "", getGeoMesh().sideSet});
    }
#endif
    return blockId;
  } else {
    std::cerr << "No cells to add." << std::endl;
    return -1;
  }
}

bool exoGeoMesh::addElemBlockProperty(const std::string &propName) {
  if (getGeoMesh().mesh->GetNumberOfCells() == 0) {
    return false;
  }
  return _elemBlockPropNames.insert(propName).second;
}

std::vector<std::string> exoGeoMesh::getElemBlockPropertyNames() const {
  std::vector<std::string> names;
  names.assign(_elemBlockPropNames.begin(), _elemBlockPropNames.end());
  return names;
}

int exoGeoMesh::getElemBlockProperty(const std::string &propName,
                                     int blockId) const {
  if (_elemBlockPropNames.find(propName) == _elemBlockPropNames.end()) {
    std::cerr << "No property found with name: " << propName << std::endl;
    exit(1);
  }
  auto blockIt = _elemBlocks.find(blockId);
  if (blockIt == _elemBlocks.end()) {
    std::cerr << "No block found with id: " << blockId << std::endl;
    exit(1);
  } else {
    auto it = blockIt->second.properties.find(propName);
    if (it == blockIt->second.properties.end()) {
      return 0;
    } else {
      return it->second;
    }
  }
}

void exoGeoMesh::setElemBlockProperty(const std::string &propName, int blockId,
                                      int value) {
  if (_elemBlockPropNames.find(propName) == _elemBlockPropNames.end()) {
    std::cerr << "No property found with name: " << propName << std::endl;
    exit(1);
  }
  auto blockIt = _elemBlocks.find(blockId);
  if (blockIt == _elemBlocks.end()) {
    std::cerr << "No block found with id: " << blockId << std::endl;
    exit(1);
  } else {
    blockIt->second.properties[propName] = value;
  }
}

int exoGeoMesh::getNumSideSets() const { return _sideSetNames.size(); }

const std::map<int, std::string> & exoGeoMesh::getSideSetNames() const {
  return _sideSetNames;
}

const std::string &exoGeoMesh::getSideSetName(int id) const {
  auto it = _sideSetNames.find(id);
  if (it == _sideSetNames.end()) {
    std::cerr << "No side set found with id: " << id << std::endl;
    exit(1);
  } else {
    return it->second;
  }
}

void exoGeoMesh::setSideSetName(int id, const std::string &name) {
  auto it = _sideSetNames.find(id);
  if (it == _sideSetNames.end()) {
    std::cerr << "No side set found with id: " << id << std::endl;
    exit(1);
  } else {
    it->second = name;
  }
}

std::vector<int> exoGeoMesh::getSideSetIds() const {
  std::vector<int> ids;
  ids.reserve(_sideSetNames.size());
  for (const auto &sideSet : _sideSetNames) {
    ids.emplace_back(sideSet.first);
  }
  return ids;
}

std::pair<std::vector<vtkIdType>, std::vector<int>> exoGeoMesh::getSideSet(
    int sideSetId) const {
  std::vector<vtkIdType> cellIdVec;
  std::vector<int> sideVec;
  if (_sideSetNames.find(sideSetId) == _sideSetNames.end()) {
    std::cerr << "No side set found with id: " << sideSetId << std::endl;
  }
  auto geoMesh = getGeoMesh();
  auto sideSet = geoMesh.sideSet;
  auto sideSetExoId = vtkIntArray::FastDownCast(
      sideSet.sides->GetCellData()->GetAbstractArray(getExoIdArrName()));
  auto sideSetOrigCellId = geoMesh.sideSet.getOrigCellArr();
  auto sideSetCellFaceId = geoMesh.sideSet.getCellFaceArr();
  for (vtkIdType i = 0; i < sideSet.sides->GetNumberOfCells(); ++i) {
    if (sideSetExoId->GetValue(i) == sideSetId) {
      cellIdVec.emplace_back(sideSetOrigCellId->GetTypedComponent(i, 0));
      sideVec.emplace_back(sideSetCellFaceId->GetTypedComponent(i, 0));
    }
  }
  return {cellIdVec, sideVec};
}

int exoGeoMesh::addSideSet(const std::vector<vtkIdType> &elements,
                           const std::vector<int> &sides,
                           const std::string &name) {
  if (!elements.empty() && elements.size() == sides.size()) {
    auto newId = nemAux::leastUnusedKey(_sideSetNames);
    auto gm = getGeoMesh();
    if (!gm.sideSet.sides) {
      auto sideSetPD = vtkSmartPointer<vtkPolyData>::New();
      sideSetPD->Allocate();
      sideSetPD->SetPoints(getGeoMesh().mesh->GetPoints());
      auto sideSetExoId = vtkSmartPointer<vtkIntArray>::New();
      sideSetExoId->SetName(getExoIdArrName());
      sideSetPD->GetCellData()->AddArray(sideSetExoId);
      auto sideSetEntities = vtkSmartPointer<vtkIntArray>::New();
      auto sideSetOrigCellId = vtkSmartPointer<vtkIdTypeArray>::New();
      sideSetOrigCellId->SetNumberOfComponents(2);
      auto sideSetCellFaceId = vtkSmartPointer<vtkIntArray>::New();
      sideSetCellFaceId->SetNumberOfComponents(2);
      setGeoMesh(
          {getGeoMesh().mesh,
           getGeoMesh().geo,
           getGeoMesh().link,
           {sideSetPD, sideSetEntities, sideSetOrigCellId, sideSetCellFaceId}});
    }
    auto sideSet = getGeoMesh().sideSet;
    auto sideSetExoId = vtkIntArray::FastDownCast(
        sideSet.sides->GetCellData()->GetAbstractArray(getExoIdArrName()));
    sideSetExoId->Resize(sideSet.sides->GetNumberOfCells() + elements.size());
    auto sideSetEntities = sideSet.getGeoEntArr();
    sideSetEntities->Resize(sideSet.sides->GetNumberOfCells() +
                            elements.size());
    auto sideSetOrigCellId = sideSet.getOrigCellArr();
    sideSetOrigCellId->Resize(sideSet.sides->GetNumberOfCells() +
                              elements.size());
    auto sideSetCellFaceId = sideSet.getCellFaceArr();
    sideSetCellFaceId->Resize(sideSet.sides->GetNumberOfCells() +
                              elements.size());

    auto mesh = getGeoMesh().mesh;
    auto elemIt = elements.begin();
    auto sidesIt = sides.begin();
    for (; elemIt != elements.end(); ++sidesIt, ++elemIt) {
      auto parent = mesh->GetCell(*elemIt);
      if (!parent) {
        std::cerr << "No cell found with index " << *elemIt << std::endl;
        continue;
      }
      auto side = parent->GetCellDimension() == 2 ? parent->GetEdge(*sidesIt)
                                                  : parent->GetFace(*sidesIt);
      sideSet.sides->InsertNextCell(side->GetCellType(), side->GetPointIds());
      sideSetExoId->InsertNextValue(newId);
      sideSetEntities->InsertNextValue(newId);
      sideSetOrigCellId->InsertNextTuple2(*elemIt, -1);
      sideSetCellFaceId->InsertNextTuple2(*sidesIt, -1);
    }

    _sideSetNames[newId] = name;
    return newId;
  }
  return -1;
}

int exoGeoMesh::getNumNodeSets() const { return _nodeSets.size(); }

std::map<int, std::string> exoGeoMesh::getNodeSetNames() const {
  std::map<int, std::string> names;
  std::transform(_nodeSets.begin(), _nodeSets.end(),
                 std::inserter(names, names.begin()),
                 [](const std::pair<int, nodeSet> &set) {
                   return std::make_pair(set.first, set.second.name);
                 });
  return names;
}

const std::string &exoGeoMesh::getNodeSetName(int id) const {
  auto it = _nodeSets.find(id);
  if (it == _nodeSets.end()) {
    std::cerr << "No side set found with id: " << id << std::endl;
    exit(1);
  } else {
    return it->second.name;
  }
}

void exoGeoMesh::setNodeSetName(int id, const std::string &name) {
  auto it = _nodeSets.find(id);
  if (it == _nodeSets.end()) {
    std::cerr << "No side set found with id: " << id << std::endl;
    exit(1);
  } else {
    it->second.name = name;
  }
}

std::vector<int> exoGeoMesh::getNodeSetIds() const {
  std::vector<int> ids;
  ids.reserve(_nodeSets.size());
  for (const auto &nodeSet : _nodeSets) {
    ids.emplace_back(nodeSet.first);
  }
  return ids;
}

const std::vector<vtkIdType> &exoGeoMesh::getNodeSet(int nodeSetId) const {
  auto it = _nodeSets.find(nodeSetId);
  if (it == _nodeSets.end()) {
    std::cerr << "No node set found with id: " << nodeSetId << std::endl;
    exit(1);
  }
  return it->second.nodes;
}

int exoGeoMesh::addNodeSet(const std::vector<vtkIdType> &nodes,
                           const std::string &name) {
  if (!nodes.empty()) {
    auto minmax = std::minmax_element(nodes.begin(), nodes.end());
    if (*minmax.first < 0) {
      std::cerr << "Unknown node " << *minmax.first << std::endl;
      return -1;
    } else if (*minmax.second >= getGeoMesh().mesh->GetNumberOfPoints()) {
      std::cerr << "Unknown node " << *minmax.second << std::endl;
      return -1;
    }
    auto newId = nemAux::leastUnusedKey(_nodeSets);
    _nodeSets[newId] = {name, nodes};
    return newId;
  }
  return -1;
}

const std::string &exoGeoMesh::getPhysGrpPropertyName() const {
  return _physGrpName;
}

void exoGeoMesh::setPhysGrpPropertyName(const std::string &physGrpName) {
  if (_elemBlockPropNames.find(physGrpName) == _elemBlockPropNames.end()) {
    std::cerr << "No property found with name " << physGrpName << std::endl;
    return;
  }
  _physGrpName = physGrpName;
  auto link = getGeoMesh().link;
  if (!link.empty()) {
    auto mesh = getGeoMesh().mesh;
    auto elemBlockIdArr = vtkIntArray::FastDownCast(
        mesh->GetCellData()->GetAbstractArray(getExoIdArrName()));
    auto linkArr = mesh->GetCellData()->GetArray(link.c_str());
    if (linkArr) {
      for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
        linkArr->SetComponent(
            i, 0,
            getElemBlockProperty(_physGrpName, elemBlockIdArr->GetValue(i)));
      }
    }
  }
}

void exoGeoMesh::stitch(const exoGeoMesh &otherGM, float tol) {
  auto mesh = this->getGeoMesh().mesh;
  if (mesh->GetNumberOfCells() == 0) {
    std::cerr << "Current mesh is empty. Cannot stitch onto empty mesh."
              << std::endl;
    exit(1);
  }
  if (this->getGeoMesh().getDimension() !=
      otherGM.getGeoMesh().getDimension()) {
    std::cerr << "Dimension of meshes does not match." << std::endl;
    exit(1);
  }
  auto numOrigCells = mesh->GetNumberOfCells();
  auto otherMesh = otherGM.getGeoMesh().mesh;
  std::map<int, int> old2newElemBlock;

  std::set<std::string> propsToCopy;
  std::set_intersection(_elemBlockPropNames.begin(), _elemBlockPropNames.end(),
                        otherGM._elemBlockPropNames.begin(),
                        otherGM._elemBlockPropNames.end(),
                        std::inserter(propsToCopy, propsToCopy.begin()));
  auto newBlockId = nemAux::leastUnusedKey(this->_elemBlocks);
  for (const auto &otherElemBlock : otherGM._elemBlocks) {
    old2newElemBlock[otherElemBlock.first] = newBlockId;
    auto &newElemBlock = this->_elemBlocks[newBlockId];
    newElemBlock.name = otherElemBlock.second.name;
    newElemBlock.cellType = otherElemBlock.second.cellType;
    for (const auto &property : propsToCopy) {
      auto it = otherElemBlock.second.properties.find(property);
      if (it != otherElemBlock.second.properties.end()) {
        newElemBlock.properties[property] = it->second;
      }
    }
    newBlockId = nemAux::leastUnusedKey(this->_elemBlocks, newBlockId);
  }
  auto oldOtherBlockIdArr =
      vtkSmartPointer<vtkIntArray>(vtkIntArray::FastDownCast(
          otherMesh->GetCellData()->GetAbstractArray(getExoIdArrName())));
  otherMesh->GetCellData()->RemoveArray(getExoIdArrName());
  auto newOtherBlockIdArr = vtkSmartPointer<vtkIntArray>::Take(
      vtkIntArray::FastDownCast(mesh->GetCellData()
                                    ->GetAbstractArray(getExoIdArrName())
                                    ->NewInstance()));
  newOtherBlockIdArr->SetName(getExoIdArrName());
  newOtherBlockIdArr->SetNumberOfValues(otherMesh->GetNumberOfCells());
  for (vtkIdType i = 0; i < otherMesh->GetNumberOfCells(); ++i) {
    newOtherBlockIdArr->SetComponent(
        i, 0, old2newElemBlock.at(oldOtherBlockIdArr->GetValue(i)));
  }
  otherMesh->GetCellData()->AddArray(newOtherBlockIdArr);

#ifdef HAVE_GMSH
  if (!this->getGeoMesh().geo.empty()) {
    gmsh::model::setCurrent(this->getGeoMesh().geo);
    gmsh::model::remove();
    mesh->GetCellData()->RemoveArray(this->getGeoMesh().link.c_str());
    this->setGeoMesh({mesh, "", "", this->getGeoMesh().sideSet});
  }
#endif
  auto newMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  auto nodeMaps = mergeMeshes(mesh, otherMesh, newMesh, tol);
  resetNodeSetPoints(nodeMaps.first);
  if (!otherGM._nodeSets.empty()) {
    auto nextNodeSetId = nemAux::leastUnusedKey(_nodeSets);
    for (const auto &otherNodeSet : otherGM._nodeSets) {
      auto &nodeSet = _nodeSets[nextNodeSetId];
      nodeSet.name = otherNodeSet.second.name;
      nodeSet.nodes.reserve(otherNodeSet.second.nodes.size());
      std::transform(
          otherNodeSet.second.nodes.begin(), otherNodeSet.second.nodes.end(),
          std::back_inserter(nodeSet.nodes),
          [&nodeMaps](vtkIdType x) { return nodeMaps.second->GetValue(x); });
      nextNodeSetId = nemAux::leastUnusedKey(_nodeSets, nextNodeSetId);
    }
  }
  resetSideSetPoints(newMesh, getGeoMesh().sideSet.sides, nodeMaps.first);
  setGeoMesh(
      {newMesh, getGeoMesh().geo, getGeoMesh().link, getGeoMesh().sideSet});
  auto otherGeoMesh = otherGM.getGeoMesh();
  if (otherGeoMesh.sideSet.sides) {
    auto otherSideSet = otherGeoMesh.sideSet;
    auto otherExoIdArr = vtkIntArray::FastDownCast(
        otherSideSet.sides->GetCellData()->GetAbstractArray(getExoIdArrName()));
    auto otherEntitiesArr = otherGeoMesh.sideSet.getGeoEntArr();
    auto otherOrigCellIdArr = otherGeoMesh.sideSet.getOrigCellArr();
    auto otherCellFaceIdArr = otherGeoMesh.sideSet.getCellFaceArr();
    auto geoMesh = this->getGeoMesh();
    auto sideSet = geoMesh.sideSet;
    vtkSmartPointer<vtkIntArray> sideSetExoId;
    vtkSmartPointer<vtkIntArray> sideSetEntities;
    vtkSmartPointer<vtkIdTypeArray> sideSetOrigCellId;
    vtkSmartPointer<vtkIntArray> sideSetCellFaceId;
    if (sideSet.sides) {
      sideSetExoId = vtkIntArray::FastDownCast(
          sideSet.sides->GetCellData()->GetAbstractArray(getExoIdArrName()));
      sideSetEntities = geoMesh.sideSet.getGeoEntArr();
      sideSetOrigCellId = geoMesh.sideSet.getOrigCellArr();
      sideSetCellFaceId = geoMesh.sideSet.getCellFaceArr();
    } else {
      auto sideSetPD = vtkSmartPointer<vtkPolyData>::New();
      sideSetPD->Allocate(otherSideSet.sides);
      sideSetExoId = vtkSmartPointer<vtkIntArray>::New();
      sideSetExoId->SetName(getExoIdArrName());
      sideSetExoId->Allocate(otherSideSet.sides->GetNumberOfCells());
      sideSet.sides->GetCellData()->AddArray(sideSetExoId);
      sideSetEntities = vtkSmartPointer<vtkIntArray>::New();
      sideSetEntities->Allocate(otherSideSet.sides->GetNumberOfCells());
      sideSetOrigCellId = vtkSmartPointer<vtkIdTypeArray>::New();
      sideSetOrigCellId->SetNumberOfComponents(2);
      sideSetOrigCellId->Allocate(otherSideSet.sides->GetNumberOfCells());
      sideSetCellFaceId = vtkSmartPointer<vtkIntArray>::New();
      sideSetCellFaceId->SetNumberOfComponents(2);
      sideSetCellFaceId->Allocate(otherSideSet.sides->GetNumberOfCells());
      sideSet = {sideSetPD, sideSetEntities, sideSetOrigCellId,
                 sideSetCellFaceId};
      this->setGeoMesh({newMesh, this->getGeoMesh().geo,
                        this->getGeoMesh().link, sideSet});
    }
    std::map<int, int> old2newSideSet;
    auto newSideSetId = nemAux::leastUnusedKey(this->_sideSetNames);
    for (const auto &otherSet : otherGM._sideSetNames) {
      old2newSideSet[otherSet.first] = newSideSetId;
      this->_sideSetNames[newSideSetId] = otherSet.second;
      newSideSetId = nemAux::leastUnusedKey(this->_sideSetNames, newSideSetId);
    }
    std::map<int, int> old2newGeoEnt;
    for (vtkIdType i = 0; i < sideSet.sides->GetNumberOfCells(); ++i) {
      old2newGeoEnt[sideSetEntities->GetValue(i)] = -1;
    }
    int nextSideSetEntId = nemAux::leastUnusedKey(old2newGeoEnt);
    for (vtkIdType i = 0; i < otherSideSet.sides->GetNumberOfCells(); ++i) {
      auto origCellId = otherOrigCellIdArr->GetTypedComponent(i, 0) + numOrigCells;
      auto cellFaceId = otherCellFaceIdArr->GetTypedComponent(i, 0);
      auto parent = newMesh->GetCell(origCellId);
      auto side = parent->GetCellDimension() == 2 ? parent->GetEdge(cellFaceId)
                                                  : parent->GetFace(cellFaceId);
      sideSet.sides->InsertNextCell(side->GetCellType(), side->GetPointIds());
      int newEntId = otherEntitiesArr->GetValue(i);
      auto insertIter = old2newGeoEnt.insert({newEntId, nextSideSetEntId});
      if (insertIter.second || insertIter.first->second == -1) {
        if (insertIter.first->second == -1) {
          insertIter.first->second = nextSideSetEntId;
        }
        nextSideSetEntId =
            nemAux::leastUnusedKey(old2newGeoEnt, nextSideSetEntId);
      }
      sideSetExoId->InsertNextValue(old2newSideSet[otherExoIdArr->GetValue(i)]);
      sideSetEntities->InsertNextValue(insertIter.first->second);
      sideSetOrigCellId->InsertNextTuple2(origCellId, -1);
      sideSetCellFaceId->InsertNextTuple2(cellFaceId, -1);
    }
  }
  otherMesh->GetCellData()->RemoveArray(getExoIdArrName());
  otherMesh->GetCellData()->AddArray(oldOtherBlockIdArr);
}

void exoGeoMesh::scaleNodes(double scale) {
  auto transform = vtkSmartPointer<vtkTransform>::New();
  transform->Scale(scale, scale, scale);
  auto transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
  transformFilter->SetTransform(transform);
  transformFilter->SetInputData(getGeoMesh().mesh);
  transformFilter->Update();
  auto mesh = vtkUnstructuredGrid::SafeDownCast(transformFilter->GetOutput());
  auto sideSet = getGeoMesh().sideSet;
  if (sideSet.sides) {
    sideSet.sides->SetPoints(mesh->GetPoints());
  }
  setGeoMesh({mesh, getGeoMesh().geo, getGeoMesh().link, getGeoMesh().sideSet});
}

geoMeshBase::GeoMesh exoGeoMesh::exoReader2GM(
    vtkSmartPointer<vtkExodusIIReader> reader) {
  auto mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  SideSet sideSet{};
  if (!reader) {
    return {mesh, "", "", {}};
  }

  auto mbds = reader->GetOutput();
  if (mbds->GetNumberOfBlocks() != 8) {
    std::cerr << "Malformed output of vtkExodusIIReader. Check that file "
                 "exists and is valid Exodus II file."
              << std::endl;
    exit(1);
  }
  // Map from "GlobalElementId", set by vtkExodusIIReader, to the index in
  // mesh, which becomes GeoMesh.mesh
  std::map<vtkIdType, vtkIdType> elemId2MeshIdx;
  if (reader->GetNumberOfElementsInFile() > 0) {
    auto exoElemBlocks = vtkMultiBlockDataSet::SafeDownCast(mbds->GetBlock(0));
    auto merge = vtkSmartPointer<mergeCells>::New();
    merge->SetUnstructuredGrid(mesh);
    merge->SetTotalNumberOfDataSets(exoElemBlocks->GetNumberOfBlocks());
    merge->SetTotalNumberOfPoints(reader->GetNumberOfNodes());
    merge->SetTotalNumberOfCells(reader->GetTotalNumberOfElements());
    merge->SetUseGlobalIds(1);
    merge->SetUseGlobalCellIds(1);
    auto it = vtkSmartPointer<vtkCompositeDataIterator>::Take(
        exoElemBlocks->NewIterator());
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextItem()) {
      merge->MergeDataSet(vtkDataSet::SafeDownCast(it->GetCurrentDataObject()));
    }
    merge->Finish();
    for (int i = mesh->GetCellData()->GetNumberOfArrays() - 1; i >= 0; --i) {
      if (mesh->GetCellData()->GetArray(i)->GetNumberOfTuples() !=
          mesh->GetNumberOfCells()) {
        mesh->GetCellData()->RemoveArray(i);
      }
    }
    for (int i = mesh->GetPointData()->GetNumberOfArrays() - 1; i >= 0; --i) {
      if (mesh->GetPointData()->GetArray(i)->GetNumberOfTuples() !=
          mesh->GetNumberOfPoints()) {
        mesh->GetPointData()->RemoveArray(i);
      }
    }
    auto globalElemIdArr = vtkIdTypeArray::FastDownCast(
        mesh->GetCellData()->GetArray(reader->GetGlobalElementIdArrayName()));
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
      elemId2MeshIdx[globalElemIdArr->GetValue(i)] = i;
    }
    mesh->GetCellData()->RemoveArray(reader->GetPedigreeElementIdArrayName());
    mesh->GetCellData()->RemoveArray(reader->GetGlobalElementIdArrayName());
  }
  if (reader->GetNumberOfSideSetArrays() > 0 &&
      reader->GetDimensionality() >= 2) {
    auto sideSetPD = vtkSmartPointer<vtkPolyData>::New();
    sideSetPD->SetPoints(mesh->GetPoints());
    sideSetPD->Allocate();
    auto sideSetExoId = vtkSmartPointer<vtkIntArray>::New();
    sideSetExoId->SetName(getExoIdArrName());
    auto sideSetEntities = vtkSmartPointer<vtkIntArray>::New();
    auto sideSetOrigCellId = vtkSmartPointer<vtkIdTypeArray>::New();
    sideSetOrigCellId->SetNumberOfComponents(2);
    auto sideSetCellFaceId = vtkSmartPointer<vtkIntArray>::New();
    sideSetCellFaceId->SetNumberOfComponents(2);

    int comp_ws = sizeof(double);
    int io_ws = 0;
    float version;
    auto exoid =
        ex_open(reader->GetFileName(), EX_READ, &comp_ws, &io_ws, &version);
    checkExodusErr(exoid, true);

    auto exoSideSets = vtkMultiBlockDataSet::SafeDownCast(mbds->GetBlock(4));
    for (unsigned int i = 0; i < exoSideSets->GetNumberOfBlocks(); ++i) {
      auto exoSideSet =
          vtkUnstructuredGrid::SafeDownCast(exoSideSets->GetBlock(i));
      // vtkExodusIIReader sets the element id off by 1 compared to
      // GlobalElementId (see vtkExodusIIReader.cxx::832 in v7.1.0)
      auto sourceElemArr = vtkIdTypeArray::FastDownCast(
          exoSideSet->GetCellData()->GetAbstractArray(
              reader->GetSideSetSourceElementIdArrayName()));
      auto objectIdArr = exoSideSet->GetCellData()->GetArray(getExoIdArrName());
      // Because the element ids are off by 1 compared to GlobalElementId,
      // the side number is off because the wrong parent cell type is used to
      // convert Exodus side index to VTK side index. Thus we do not use
      // exoSideSet->GetCellData()->GetAbstractArray(
      //     reader->GetSideSetSourceElementSideArrayName());
      // and instead read from the file.
      // Here we rely on the assumption that vtkExodusIIReader does not change
      // the order of entries in the side set.
      std::vector<int> sourceSideVec;
      sourceSideVec.resize(exoSideSet->GetNumberOfCells());
      auto err = ex_get_set(
          exoid, EX_SIDE_SET,
          vtkIntArray::FastDownCast(
              exoSideSet->GetCellData()->GetAbstractArray(getExoIdArrName()))
              ->GetValue(0),
          nullptr, sourceSideVec.data());
      checkExodusErr(err);

      for (vtkIdType j = 0; j < exoSideSet->GetNumberOfCells(); ++j) {
        auto origCellId = elemId2MeshIdx.at(sourceElemArr->GetValue(j) + 1);
        auto parent = mesh->GetCell(origCellId);
        auto cellFaceId = exoSide2vtkSide(
            sourceSideVec[j], static_cast<VTKCellType>(parent->GetCellType()));
        if (reader->GetDimensionality() == parent->GetCellDimension()) {
          if (reader->GetDimensionality() == 2) {
            auto edge = parent->GetEdge(cellFaceId);
            sideSetPD->InsertNextCell(edge->GetCellType(), edge->GetPointIds());
          } else if (reader->GetDimensionality() == 3) {
            auto face = parent->GetFace(cellFaceId);
            sideSetPD->InsertNextCell(face->GetCellType(), face->GetPointIds());
          }
          sideSetExoId->InsertNextValue(
              static_cast<int>(objectIdArr->GetComponent(j, 0)));
          sideSetEntities->InsertNextValue(
              static_cast<int>(objectIdArr->GetComponent(j, 0)));
          sideSetOrigCellId->InsertNextTuple2(origCellId, -1);
          sideSetCellFaceId->InsertNextTuple2(cellFaceId, -1);
        }
      }
    }
    sideSetPD->GetCellData()->AddArray(sideSetExoId);
    sideSetPD->Squeeze();
    sideSet = {sideSetPD, sideSetEntities, sideSetOrigCellId,
               sideSetCellFaceId};
    auto err = ex_close(exoid);
    checkExodusErr(err);
  }
  return {mesh, "", "", sideSet};
}

vtkSmartPointer<vtkExodusIIReader> exoGeoMesh::getExoReader(
    const std::string &fileName, int timeStep) {
  static constexpr std::array<vtkExodusIIReader::ObjectType, 8> objectTypes = {
      vtkExodusIIReader::ObjectType::ELEM_BLOCK_ELEM_CONN,
      vtkExodusIIReader::ObjectType::NODE_SET_CONN,
      vtkExodusIIReader::ObjectType::SIDE_SET_CONN,
      vtkExodusIIReader::ObjectType::NODAL,
      vtkExodusIIReader::ObjectType::GLOBAL,
      vtkExodusIIReader::ObjectType::ELEM_BLOCK,
      vtkExodusIIReader::ObjectType::NODE_SET,
      vtkExodusIIReader::ObjectType::SIDE_SET};
  if (fileName.empty()) {
    return nullptr;
  }
  auto reader = vtkSmartPointer<vtkExodusIIReader>::New();
  reader->SetFileName(fileName.c_str());
  // time steps are 1-indexed in Exodus
  reader->SetTimeStep(timeStep - 1);
  reader->UpdateInformation();
  for (const auto &objectType : objectTypes) {
    reader->SetAllArrayStatus(objectType, 1);
  }
  for (int i = 0; i < reader->GetNumberOfNodeMapArrays(); ++i) {
    reader->SetNodeMapArrayStatus(reader->GetNodeMapArrayName(i), 0);
  }
  for (int i = 0; i < reader->GetNumberOfElementMapArrays(); ++i) {
    reader->SetElementMapArrayStatus(reader->GetElementMapArrayName(i), 0);
  }
  reader->GenerateGlobalElementIdArrayOn();
  reader->GenerateGlobalNodeIdArrayOn();
  reader->Update();
  return reader;
}

void exoGeoMesh::resetNative() {
  _nodeSets.clear();
  auto mesh = this->getGeoMesh().mesh;
  auto link = this->getGeoMesh().link;
  _physGrpName = link;
  if (!mesh->GetCellData()->HasArray(getExoIdArrName())) {
    auto elemBlockIdArr = vtkSmartPointer<vtkIntArray>::New();
    elemBlockIdArr->SetName(getExoIdArrName());
    auto linkArr = mesh->GetCellData()->GetArray(link.c_str());
    if (!link.empty() && linkArr) {
      elemBlockIdArr->Allocate(mesh->GetNumberOfCells());
      for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); ++i) {
        elemBlockIdArr->InsertNextValue(
            static_cast<int>(linkArr->GetComponent(i, 0)));
      }
    } else {
      elemBlockIdArr->SetNumberOfValues(mesh->GetNumberOfCells());
      elemBlockIdArr->FillComponent(0, 1);
    }
    mesh->GetCellData()->AddArray(elemBlockIdArr);
  }
  resetElemBlocks();
  auto gm = getGeoMesh();
  gm.findSide2OrigCell();
  setGeoMesh(gm);
  setSideSetObjId();
}

void exoGeoMesh::resetElemBlocks() {
  _elemBlocks.clear();
  _elemBlockPropNames.clear();
  if (!_physGrpName.empty()) {
    _elemBlockPropNames.insert(_physGrpName);
  }
  std::map<std::pair<int, VTKCellType>, int> old2newElemBlock;
  auto mesh = getGeoMesh().mesh;
  auto elemBlockIds = vtkIntArray::SafeDownCast(
      mesh->GetCellData()->GetAbstractArray(getExoIdArrName()));
  vtkIdType i = 0;
  auto nextElemBlockId = nemAux::leastUnusedKey(_elemBlocks);
  for (auto it =
           vtkSmartPointer<vtkCellIterator>::Take(mesh->NewCellIterator());
       !it->IsDoneWithTraversal(); it->GoToNextCell()) {
    auto oldElemBlockId = elemBlockIds->GetValue(i);
    auto cellType = static_cast<VTKCellType>(it->GetCellType());
    auto inserted =
        old2newElemBlock.insert({{oldElemBlockId, cellType}, nextElemBlockId});
    if (inserted.second) {
      auto &elemBlock = _elemBlocks[nextElemBlockId];
      elemBlock.cellType = cellType;
      if (!_physGrpName.empty()) {
        elemBlock.properties[_physGrpName] = oldElemBlockId;
      }
      nextElemBlockId = nemAux::leastUnusedKey(_elemBlocks);
    }
    elemBlockIds->SetValue(i, inserted.first->second);
    ++i;
  }
}

void exoGeoMesh::resetNodeSetPoints(vtkIdTypeArray *nodeMap) {
  for (auto &nodeSet : _nodeSets) {
    std::unordered_set<vtkIdType> set;
    for (const auto &node : nodeSet.second.nodes) {
      set.insert(nodeMap ? nodeMap->GetValue(node) : node);
    }
    nodeSet.second.nodes.assign(set.begin(), set.end());
  }
}

void exoGeoMesh::setSideSetObjId() {
  _sideSetNames.clear();
  auto gm = getGeoMesh();
  auto mesh = gm.mesh;
  auto link = gm.link;
  auto sideSet = gm.sideSet;
  if (sideSet.sides) {
    auto sideSetObjId = vtkIntArray::FastDownCast(
        sideSet.sides->GetCellData()->GetAbstractArray(getExoIdArrName()));
    if (sideSetObjId) {
      for (vtkIdType i = 0; i < sideSet.sides->GetNumberOfCells(); ++i) {
        _sideSetNames[sideSetObjId->GetValue(i)];
      }
    } else {
      sideSetObjId = vtkIntArray::New();
      sideSetObjId->SetName(getExoIdArrName());
      sideSetObjId->Allocate(sideSet.sides->GetNumberOfCells());
      sideSet.sides->GetCellData()->AddArray(sideSetObjId);
      auto sideSetEntities = gm.sideSet.getGeoEntArr();
      auto nameArr = sideSet.getSideSetNames();
      // Exodus II requires IDs to be positive, which is not always the case for
      // other mesh types, hence use of extra array
      std::map<int, int> geoEnt2ExoId;
      int nextSideSetId = 1;
      for (vtkIdType i = 0; i < sideSet.sides->GetNumberOfCells(); ++i) {
        int geoEnt = sideSetEntities->GetValue(i);
        // If geoEnt is not found in geoEnt2ExoID, call to operator[] to sets
        // geoEnt2ExoId[geoEnt] to 0 (which is used as sentinel value)
        auto &exoID = geoEnt2ExoId[geoEnt];
        if (exoID == 0) {
          exoID = nextSideSetId;
          auto emplaceIter = _sideSetNames.emplace(exoID, std::string{});
          if (emplaceIter.second && nameArr && geoEnt >= 0 &&
              geoEnt < nameArr->GetNumberOfValues()) {
            emplaceIter.first->second = nameArr->GetValue(geoEnt);
          }
          ++nextSideSetId;
        }
        sideSetObjId->InsertNextValue(exoID);
      }
      sideSetObjId->Delete();
    }
  }
}

}  // namespace MSH
}  // namespace NEM
