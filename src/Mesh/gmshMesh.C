#include <gmsh.h>
#include <gmsh/Context.h>
#include <gmsh/CreateFile.h>
#include <gmsh/MHexahedron.h>
#include <gmsh/MLine.h>
#include <gmsh/MPoint.h>
#include <gmsh/MPrism.h>
#include <gmsh/MPyramid.h>
#include <gmsh/MQuadrangle.h>
#include <gmsh/MTetrahedron.h>
#include <gmsh/MTriangle.h>

#include "Mesh/gmshMesh.H"

#include <vtkIdList.h>
#include <vtkUnstructuredGrid.h>

class nemosysGModel : public GModel {
 public:
  void storeElementsInEntities(std::map<int, std::vector<MElement *>> &map) {
    GModel::_storeElementsInEntities(map);
  }

  void associateEntityWithMeshVertices(bool force = false) {
    GModel::_associateEntityWithMeshVertices(force);
  }

  void storeVerticesInEntities(std::vector<MVertex *> &vertices) {
    GModel::_storeVerticesInEntities(vertices);
  }

  void storePhysicalTagsInEntities(
      int dim, std::map<int, std::map<int, std::string>> &map) {
    GModel::_storePhysicalTagsInEntities(dim, map);
  }
};

gmshMesh::gmshMesh() {
  gmsh::initialize();
  //  _gmshGModel = GModel();
  std::cout << "gmshMesh constructed" << std::endl;
}

gmshMesh::gmshMesh(meshBase *mb) {
  // this constructor tries a brute-force
  // method. Use with care.
  gmsh::initialize();
  // dataSet uses smart pointers so deep
  // copy is not needed
  dataSet = mb->getDataSet();
  numPoints = dataSet->GetNumberOfPoints();
  numCells = dataSet->GetNumberOfCells();
  std::cout << "Number of points: " << numPoints << std::endl;
  std::cout << "Number of cells: " << numCells << std::endl;
  std::cout << "gmshMesh constructed" << std::endl;
}

gmshMesh::~gmshMesh() {
  //  _gmshGModel.destroy();
  gmsh::finalize();
  std::cout << "gmshMesh destroyed" << std::endl;
}

gmshMesh::gmshMesh(const std::string &fname) {
  gmsh::initialize();
  filename = fname;
  read(filename);
  std::cout << "gmshMesh constructed" << std::endl;
}

void gmshMesh::read(const std::string &fname) {
  // Read the file.
  GModel _gmshGModel = GModel();
  _gmshGModel.load(fname);

  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  // get all entities in the model
  std::vector<GEntity *> entities;
  _gmshGModel.getEntities(entities);

  // Add mesh vertices
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  // Map from gmsh Index `long int` to vtk Id `vtkIdType`. Necessary since gmsh
  // indices are not their position in the `mesh_vertices` array.
  std::map<long int, vtkIdType> gmshIndex2vtkId;
  for (const auto &entity : entities)
    for (const auto &mesh_vertex : entity->mesh_vertices) {
      // Order of inserting into map and inserting into vtkPoints is important.
      // Be mindful if editing.
      gmshIndex2vtkId[mesh_vertex->getIndex()] = points->GetNumberOfPoints();
      // `SPoint3` has implicit conversion into `double *`
      points->InsertNextPoint(mesh_vertex->point());

      //      std::cout << "Inserted point. From mesh_vertex: "
      //                << mesh_vertex->point().x() << "  "
      //                << mesh_vertex->point().y() << "  "
      //                << mesh_vertex->point().z() << "  "
      //                << "  To vtkPoints: "
      //                <<
      //                *(points->GetPoint(gmshIndex2vtkId[mesh_vertex->getIndex()]))
      //                << "  "
      //                <<
      //                *(points->GetPoint(gmshIndex2vtkId[mesh_vertex->getIndex()])+1)
      //                << "  "
      //                <<
      //                *(points->GetPoint(gmshIndex2vtkId[mesh_vertex->getIndex()])+2)
      //                << "  "
      //                << std::endl;
    }
  dataSet_tmp->SetPoints(points);

  // Count mesh elements
  std::size_t numMeshElements = 0;
  for (const auto &entity : entities)
    numMeshElements += entity->getNumMeshElements();
  dataSet_tmp->Allocate(numMeshElements);

  // Add mesh elements
  const char **typeName = nullptr;
  std::vector<const char **> unsupportedTypesList{};
  std::size_t numUnsupportedElements = 0;
  for (const auto &entity : entities) {
    for (int j = 0; j != entity->getNumMeshElements(); ++j) {
      // Check if element has an equivalent VTK cell type.
      if (entity->getMeshElement(j)->getTypeForVTK()) {
        vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
        // Allocate vertices
        idList->Allocate(entity->getMeshElement(j)->getNumVertices());
        // Add vertices to vtkIdList
        for (std::size_t k = 0;
             k != entity->getMeshElement(j)->getNumVertices(); ++k)
          idList->InsertNextId(
              //              entity->getMeshElement(j)->getVertexVTK(k)->getIndex()
              //              - 1);
              gmshIndex2vtkId
                  [entity->getMeshElement(j)->getVertex(k)->getIndex()]);
        // Add cell to vtkUnstructuredGrid
        dataSet_tmp->InsertNextCell(entity->getMeshElement(j)->getTypeForVTK(),
                                    idList);
      } else {
        ++numUnsupportedElements;
        MElement::getInfoMSH(entity->getMeshElement(j)->getTypeForMSH(),
                             typeName);
        unsupportedTypesList.emplace_back(typeName);
      }
    }
  }

  if (numUnsupportedElements) {
    std::cerr << "WARNING: " << numUnsupportedElements
              << " unsupported Gmsh elements detected of the following types:";
    for (const auto &unsupportedType : unsupportedTypesList) {
      std::cerr << " \"" << unsupportedType << "\"";
    }
    std::cerr << std::endl;
  }

  /* TODO:
   * Add geometry storage as vtkCellData
   * */

  dataSet = dataSet_tmp;

  numPoints = dataSet->GetNumberOfPoints();
  numCells = dataSet->GetNumberOfCells();
  std::cout << "Number of points: " << numPoints << std::endl;
  std::cout << "Number of cells: " << numCells << std::endl;
}
void gmshMesh::write(const std::string &fname) const {
  write(fname, 4.1, false);
}

void gmshMesh::write(const std::string &fname, double mshFileVersion,
                     bool binary) const {
  nemosysGModel _nemosysGModel = nemosysGModel();

  // Add mesh vertices
  std::vector<MVertex *> vertices(numPoints);
  double xyz[3];  // Reused in mesh element loop.

  for (nemId_t i = 0; i != numPoints; ++i) {
    dataSet->GetPoint(i, xyz);
    vertices[i] = new MVertex(xyz[0], xyz[1], xyz[2]);
  }

  // Add mesh elements
  std::map<int, std::vector<MElement *>> elements[8];

  for (nemId_t i = 0; i != numCells; ++i) {
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    // Assign type
    int type = dataSet->GetCellType(i);

    dataSet->GetCellPoints(i, idList);

    std::vector<MVertex *> cell(idList->GetNumberOfIds());
    for (vtkIdType j = 0; j != idList->GetNumberOfIds(); ++j) {
      cell[j] = vertices[idList->GetId(j)];
    }

    /* vv ------- Taken from Gmsh source: src/Geo/GModelIO_VTK.cpp:278 -------
     * vv */
    switch (type) {
      case 1:
        elements[0][1].push_back(new MPoint(cell));
        break;
      // first order elements
      case 3:
        elements[1][1].push_back(new MLine(cell));
        break;
      case 5:
        elements[2][1].push_back(new MTriangle(cell));
        break;
      case 9:
        elements[3][1].push_back(new MQuadrangle(cell));
        break;
      case 10:
        elements[4][1].push_back(new MTetrahedron(cell));
        break;
      case 12:
        elements[5][1].push_back(new MHexahedron(cell));
        break;
      case 13:
        elements[6][1].push_back(new MPrism(cell));
        break;
      case 14:
        elements[7][1].push_back(new MPyramid(cell));
        break;
      // second order elements
      case 21:
        elements[1][1].push_back(new MLine3(cell));
        break;
      case 22:
        elements[2][1].push_back(new MTriangle6(cell));
        break;
      case 23:
        elements[3][1].push_back(new MQuadrangle8(cell));
        break;
      case 28:
        elements[3][1].push_back(new MQuadrangle9(cell));
        break;
      case 24:
        elements[4][1].push_back(new MTetrahedron10(cell));
        break;
      case 25:
        elements[5][1].push_back(new MHexahedron20(cell));
        break;
      case 29:
        elements[5][1].push_back(new MHexahedron27(cell));
        break;
      default:
        Msg::Error("Unknown type of cell %d", type);
        break;
    }
    /* ^^ ------- Taken from Gmsh source: src/Geo/GModelIO_VTK.cpp:278 -------
     * ^^ */
  }

  for (auto &&element : elements)
    _nemosysGModel.storeElementsInEntities(element);
  _nemosysGModel.associateEntityWithMeshVertices();
  _nemosysGModel.storeVerticesInEntities(vertices);

  // store the physical tags
  //  for(int i = 0; i < 4; i++)
  //    _nemosysGModel.storePhysicalTagsInEntities(i, physicals[i]);
  // Renumber the indices to be contiguous.
  _nemosysGModel.renumberMeshVertices();
  _nemosysGModel.renumberMeshElements();

  // Write the file.
  //  _nemosysGModel.save(fname);

  //  GModel *temp = nemosysGModel::current();
  //  _nemosysGModel.setAsCurrent();
  //  CTX::instance()->mesh.
  //  CreateOutputFile(fname, FORMAT_MSH);
  //  nemosysGModel::setCurrent(temp);
  _nemosysGModel.writeMSH(
      fname, mshFileVersion, binary  //,
      //                          saveAll, saveParametric, scalingFactor
  );
  std::cout << "gmshMesh saved to " << fname << std::endl;
}

std::vector<double> gmshMesh::getPoint(nemId_t id) const {
  return std::vector<double>();
}
std::vector<std::vector<double>> gmshMesh::getVertCrds() const {
  return std::vector<std::vector<double>>();
}
std::map<nemId_t, std::vector<double>> gmshMesh::getCell(nemId_t id) const {
  return std::map<nemId_t, std::vector<double>>();
}
std::vector<std::vector<double>> gmshMesh::getCellVec(nemId_t id) const {
  return std::vector<std::vector<double>>();
}
void gmshMesh::inspectEdges(const std::string &ofname) const {}
vtkSmartPointer<vtkDataSet> gmshMesh::extractSurface() {
  return vtkSmartPointer<vtkDataSet>();
}
std::vector<double> gmshMesh::getCellLengths() const {
  return std::vector<double>();
}
std::vector<double> gmshMesh::getCellCenter(nemId_t cellID) const {
  return std::vector<double>();
}
int gmshMesh::getCellType() const { return 0; }
std::vector<nemId_t> gmshMesh::getConnectivities() const {
  return std::vector<nemId_t>();
}

void gmshMesh::report() const { meshBase::report(); }
