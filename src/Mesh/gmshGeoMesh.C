#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "gmshGeoMesh.H"

#include <array>
#include <map>
#include <set>

#include <gmsh.h>

#include <vtkCellIterator.h>
#include <vtkCellTypes.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

#include "AuxiliaryFunctions.H"

namespace NEM {
namespace MSH {

vtkStandardNewMacro(gmshGeoMesh)

gmshGeoMesh *gmshGeoMesh::Read(const std::string &fileName) {
  GmshInterface::Initialize();

  gmsh::open(fileName);

  std::string gmshMesh;
  gmsh::model::getCurrent(gmshMesh);

  return new gmshGeoMesh(gmshMesh);
}

gmshGeoMesh::gmshGeoMesh() : gmshGeoMesh(std::string()) {}

gmshGeoMesh::gmshGeoMesh(const std::string &gmshMesh)
    : geoMeshBase(gmsh2GM(gmshMesh)), _gmshMesh(gmshMesh) {
  std::cout << "gmshGeoMesh constructed" << std::endl;
}

gmshGeoMesh::~gmshGeoMesh() {
  GmshInterface::Finalize();
  std::cout << "gmshGeoMesh destructed" << std::endl;
}

void gmshGeoMesh::write(const std::string &fileName) {
  _gmshMesh = GM2gmsh(this->getGeoMesh());

  /*
  {  // DEBUG
    std::cout << "_gmshMesh: " << _gmshMesh << std::endl;
    std::string name;
    gmsh::model::getCurrent(name);
    std::cout << name << std::endl;
  }
  */

  gmsh::write(fileName);
}

void gmshGeoMesh::report(std::ostream &out) const { geoMeshBase::report(out); }

// TODO: use numV to sort into higher order elements.
VTKCellType gmshGeoMesh::getVTKTypeFromGmshType(const std::string &gmshType) {
  std::string name = nemAux::findToStr(gmshType, " ");
  int numV;
  if(name == "Point") {
    numV = 1;
  } else {
    numV = std::stoi(nemAux::findFromStr(gmshType, " "));
  }

  if (name == "Point") {
    return VTK_VERTEX;
  } else if (name == "Line") {
    switch(numV) {
      case 3 : return VTK_QUADRATIC_EDGE;
      default : return VTK_LINE;
    }
  } else if (name == "Triangle") {
    switch(numV) {
      case 6 : return VTK_QUADRATIC_TRIANGLE;
      default : return VTK_TRIANGLE;
    }
  } else if (name == "Quadrilateral") {
    switch(numV) {
      case 8 : return VTK_QUADRATIC_QUAD;
      case 9 : return VTK_BIQUADRATIC_QUAD;
      default : return VTK_QUAD;
    }
  } else if (name == "Polygon") {
    return VTK_POLYGON;
  } else if (name == "Tetrahedron") {
    switch(numV) {
      case 10 : return VTK_QUADRATIC_TETRA;
      default : return VTK_TETRA;
    }
  } else if (name == "Hexahedron") {
    switch(numV) {
      case 20 : return VTK_QUADRATIC_HEXAHEDRON;
      case 27 : return VTK_TRIQUADRATIC_HEXAHEDRON;
      default : return VTK_HEXAHEDRON;
    }
  } else if (name == "Prism") {
    return VTK_WEDGE;
  } else if (name == "Pyramid") {
    return VTK_PYRAMID;
  } else if (name == "Trihedron" ||  // Flat Quad + 2 Tris.
             name == "Polyhedron"    // Occurs in MElementCuts
  ) {
    std::cerr << "ERROR in Gmsh to VTK element conversion: Gmsh element \""
              << name << "\" is not supported by VTK." << std::endl;
  } else {
    std::cerr << "ERROR in Gmsh to VTK element conversion: Gmsh element \""
              << name << "\" is not recognized." << std::endl;
  }
  throw;
}

geoMeshBase::GeoMesh gmshGeoMesh::gmsh2GM(const std::string &gmshMesh) {
  auto vtkMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

  if (gmshMesh.empty()) return {vtkMesh, gmshMesh, "", nullptr};

  std::string geoEntArrayName = GEO_ENT_DEFAULT_NAME;

  gmsh::model::setCurrent(gmshMesh);

  // Track GMSH node ids and VTK point ids.
  std::map<std::size_t, vtkIdType> gmshNodeNum2vtkId;
  // Track GMSH element ids and VTK cell ids
  std::map<std::size_t, vtkIdType> gmshElementNum2vtkId;

  {  // Add points
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::Take(vtkPoints::New(VTK_DOUBLE));

    // Get GMSH nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

    // Allocate
    points->Allocate(nodeTags.size());

    // Convert GMSH nodes to VTK points
    for (std::size_t i = 0; i < nodeTags.size(); ++i) {
      gmshNodeNum2vtkId[nodeTags[i]] = points->GetNumberOfPoints();
      points->InsertNextPoint(coord[i * 3 + 0], coord[i * 3 + 1],
                              coord[i * 3 + 2]);
    }

    vtkMesh->SetPoints(points);
  }

  {  // Add elements and geometric entities
    // Note: Cannot pre-allocate cells since GMSH does not provide an API to get
    // total number of elements

    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags);

    // geoEnt: Track entities to place as cell data
    vtkSmartPointer<vtkIntArray> vtkEntities =
        vtkSmartPointer<vtkIntArray>::New();
    vtkEntities->Initialize();  // Cannot allocate
    vtkEntities->SetName(geoEntArrayName.c_str());

    // Parse by entity
    for (const auto &dimTag : dimTags) {
      // Get GMSH elements
      std::vector<int> elementTypes;
      std::vector<std::vector<std::size_t>> elementTags;
      std::vector<std::vector<std::size_t>> nodeTags;
      gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags,
                                     dimTag.first, dimTag.second);

      // Convert GMSH elements to VTK cells
      for (std::size_t iety = 0; iety < elementTypes.size(); ++iety) {
        // Get element type information from GMSH
        std::string elementName;
        int dim, order, numNodes, numPrimaryNodes;
        std::vector<double> nodeCoord;
        gmsh::model::mesh::getElementProperties(elementTypes[iety], elementName,
                                                dim, order, numNodes, nodeCoord,
                                                numPrimaryNodes);

        // Convert GMSH element type to VTK cell type
        VTKCellType ct = getVTKTypeFromGmshType(elementName);
        if (dim == gmsh::model::getDimension()) {
          for (std::size_t ietg = 0; ietg < elementTags[iety].size(); ++ietg) {
            gmshElementNum2vtkId[elementTags[iety][ietg]] =
                vtkMesh->GetNumberOfCells();

            // Create list of point ids for the cell.
            vtkSmartPointer<vtkIdList> ptIds =
                vtkSmartPointer<vtkIdList>::New();
            ptIds->Allocate(numNodes);
            for (int nn = 0; nn < numNodes; ++nn) {
              ptIds->InsertNextId(
                  gmshNodeNum2vtkId[nodeTags[iety][ietg * numNodes + nn]]);
            }

            // Add cell to vtkUnstructuredGrid
            vtkMesh->InsertNextCell(ct, ptIds);

            // geoEnt: Add to geometric entity array
            vtkEntities->InsertNextValue(dimTag.second);
          }
        }
      }
    }

    // geoEnt: Add geometric entity data to VTK
    vtkMesh->GetCellData()->AddArray(vtkEntities);
  }

  {
      // Add point data
      // TODO
  }

  {  // Add cell data
     // TODO
  }

  // Since could not allocate elements/cells, we have to squeeze
  vtkMesh->Squeeze();

  // Only keep the geometry in the gmshMesh
  gmsh::model::mesh::clear();

  return {vtkMesh, gmshMesh, geoEntArrayName, nullptr};
}

std::string gmshGeoMesh::GM2gmsh(const GeoMesh &geoMesh) {
  /*
  // DEBUG
  std::cout << "geoMesh.geo: " << geoMesh.geo << std::endl;
  */

  gmsh::model::setCurrent(geoMesh.geo);

  /*
  {  // DEBUG
    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags);
    for (const auto &dimTag : dimTags) {
      std::cout << "geoEnt   dim: " << dimTag.first
                << "  tag: " << dimTag.second << std::endl;
    }

    gmsh::model::getPhysicalGroups(dimTags);
    for (const auto &dimTag : dimTags) {
      std::cout << "phyGrp   dim: " << dimTag.first
                << "  tag: " << dimTag.second << std::endl;
    }
  }
  */

  if (geoMesh.mesh->GetNumberOfPoints() > 0) {  // Add points
    std::vector<std::size_t> nodeTags;
    std::vector<double> coord;

    vtkSmartPointer<vtkPoints> dataSetPoints = geoMesh.mesh->GetPoints();
    coord.resize(3 * dataSetPoints->GetNumberOfPoints());
    nodeTags.resize(dataSetPoints->GetNumberOfPoints());
    for (vtkIdType ptId = 0; ptId < dataSetPoints->GetNumberOfPoints();
         ++ptId) {
      nodeTags[ptId] = ptId + 1;
      dataSetPoints->GetPoint(ptId, &coord[3 * ptId]);
    }

    gmsh::vectorpair dimTags;
    gmsh::model::getEntities(dimTags);
    gmsh::model::mesh::addNodes(dimTags[0].first, dimTags[0].second, nodeTags,
                                coord);
  }

  if (geoMesh.mesh->GetNumberOfCells() > 0) {  // Add elements
    /*
    // Debug
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> nodeTags;
    */

    // Sort all cells by entity and element type
    vtkSmartPointer<vtkCellIterator> it = geoMesh.mesh->NewCellIterator();
    vtkSmartPointer<vtkIntArray> geoEntArray = vtkIntArray::FastDownCast(
        geoMesh.mesh->GetCellData()->GetArray(geoMesh.link.c_str()));
    for (it->InitTraversal(); !it->IsDoneWithTraversal(); it->GoToNextCell()) {
      int geoEnt = geoEntArray->GetTypedComponent(it->GetCellId(), 0);

      vtkSmartPointer<vtkIdList> ptIds = it->GetPointIds();
      std::vector<std::size_t> nodeTag;
      for (vtkIdType pt = 0; pt < it->GetNumberOfPoints(); ++pt) {
        nodeTag.emplace_back(ptIds->GetId(pt) + 1);
      }

      gmsh::model::mesh::addElementsByType(
          geoEnt, getGmshTypeFromVTKType(it->GetCellType()),
          {static_cast<std::size_t>(it->GetCellId() + 1)}, nodeTag);

      /*
      {  // Debug
        std::cout << "cellId: " << it->GetCellId() + 1 << "  geoEnt: " << geoEnt
                  << "  nodeTags:";
        for (const auto &node : nodeTag) std::cout << " " << node;
        std::cout << std::endl;
      }
      */
    }
    it->Delete();
  }

  /*
  {  // DEBUG
    int NbNode;
    double NbNodeD;
    gmsh::option::getNumber("Mesh.NbNodes", NbNodeD);
    NbNode = static_cast<int>(NbNodeD);
    std::cout << "NbNode: " << NbNode << std::endl;
  }
  */

  gmsh::model::mesh::reclassifyNodes();

  return geoMesh.geo;
}

void gmshGeoMesh::resetNative() {
  auto gm = getGeoMesh();
  if (gm.geo.empty()) {
    gm.geo = "geoMesh_" + nemAux::getRandomString(6);
    gmsh::model::add(gm.geo);
  }
  if (gm.link.empty()) {
    gm.link = GEO_ENT_DEFAULT_NAME;
    auto linkArr = vtkSmartPointer<vtkIntArray>::New();
    linkArr->SetName(gm.link.c_str());
    linkArr->SetNumberOfTuples(getGeoMesh().mesh->GetNumberOfCells());
    linkArr->FillComponent(0, 1);
    getGeoMesh().mesh->GetCellData()->AddArray(linkArr);
  }
  setGeoMesh(gm);
  _gmshMesh = GM2gmsh(gm);
}

}  // namespace MSH
}  // namespace NEM
