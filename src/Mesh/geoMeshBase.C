#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif

#include "geoMeshBase.H"

#include <iostream>
#include <utility>

#include <gmsh.h>

namespace NEM {
namespace MSH {

GmshInterface::GmshInterface() {
  gmsh::initialize();
  gmsh::option::setNumber("General.Terminal", 1.0);  // Gmsh errors to stderr
  gmsh::option::setNumber("Mesh.SaveAll", 1);
  std::cout << "Gmsh initialized" << std::endl;
}

GmshInterface::~GmshInterface() {
  gmsh::finalize();
  std::cout << "Gmsh finalized" << std::endl;
}

std::shared_ptr<GmshInterface> GmshInterface::instance = nullptr;
int GmshInterface::count = 0;

void GmshInterface::Initialize() {
  ++count;
  if (!instance) instance = std::shared_ptr<GmshInterface>(new GmshInterface());
}

void GmshInterface::Finalize() {
  --count;
  if (count == 0) instance.reset();
}

geoMeshBase::geoMeshBase()
    : geoMeshBase({vtkSmartPointer<vtkUnstructuredGrid>::New(), "", ""}) {}

geoMeshBase::geoMeshBase(GeoMesh inGeoMesh)
    : _geoMesh(std::move(inGeoMesh)), _angleThreshold(30.0 * M_PI / 180.0) {
  GmshInterface::Initialize();
  std::cout << "geoMeshBase constructed" << std::endl;
}

geoMeshBase::~geoMeshBase() {
  GmshInterface::Finalize();
  std::cout << "geoMeshBase destructed" << std::endl;
}

// void geoMeshBase::update() {
//  /*
//  if (_stale.GEO_MESH_LINK && !_stale.GEOMETRY && !_stale.MESH) {
//    // Need to loop through all entities in the geometry and all points and
//    // cells in the mesh and set the entity ID.
//    std::vector<GEntity *> entities;
//    _GModel.getEntities(entities, 1);
//    // entities[0]->containsPoint()
//
//    for (vtkIdType ptId = 0; ptId != _mesh->GetNumberOfPoints(); ++ptId) {
//      std::array<double, 3> pt{};
//      _mesh->GetPoint(ptId, pt.data());
//    }
//  }
//  */
//}

}  // namespace MSH
}  // namespace NEM
