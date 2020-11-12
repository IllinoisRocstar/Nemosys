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
    : geoMeshBase(
          {vtkSmartPointer<vtkUnstructuredGrid>::New(), "", "", nullptr}) {}

geoMeshBase::geoMeshBase(GeoMesh inGeoMesh)
    : _geoMesh(std::move(inGeoMesh)), _angleThreshold(30.0 * M_PI / 180.0) {
  InitializeObjectBase();
  GmshInterface::Initialize();
  std::cout << "geoMeshBase constructed" << std::endl;
}

geoMeshBase::~geoMeshBase() {
  this->ReferenceCount--;
  GmshInterface::Finalize();
  std::cout << "geoMeshBase destructed" << std::endl;
}

void geoMeshBase::report(std::ostream &out) const {
  out << this->GetClassName() << ":\n";
  auto mesh = getGeoMesh().mesh;
  out << "Nodes:\t" << mesh->GetNumberOfPoints() << '\n';
  out << "Elements:\t" << mesh->GetNumberOfCells() << '\n';
}

void geoMeshBase::takeGeoMesh(geoMeshBase *otherGeoMesh) {
  _geoMesh = std::move(otherGeoMesh->_geoMesh);
  otherGeoMesh->_geoMesh = {
      vtkSmartPointer<vtkUnstructuredGrid>::New(), {}, {}, nullptr};
  otherGeoMesh->resetNative();
  resetNative();
}

// Helper functions for geoMeshBase::get{Point, Cell, Field}DataArrayCopy
namespace {

vtkSmartPointer<vtkAbstractArray> copyOrCreateCopy(
    vtkAbstractArray *orig, vtkAbstractArray *copy = nullptr) {
  if (copy) {
    copy->DeepCopy(orig);
    return nullptr;
  } else if (orig) {
    auto newCopy = vtkSmartPointer<vtkAbstractArray>::Take(orig->NewInstance());
    newCopy->DeepCopy(orig);
    return newCopy;
  } else {
    return nullptr;
  }
}

vtkSmartPointer<vtkAbstractArray> getArrCopy(
    vtkFieldData *field, const std::string &arrayName, int *arrayIdx,
    vtkAbstractArray *array = nullptr) {
  auto orig = arrayIdx ? field->GetAbstractArray(arrayName.c_str(), *arrayIdx)
                       : field->GetAbstractArray(arrayName.c_str());
  return copyOrCreateCopy(orig, array);
}

vtkSmartPointer<vtkAbstractArray> getArrCopy(
    vtkFieldData *field, int arrayIdx, vtkAbstractArray *array = nullptr) {
  auto orig = field->GetAbstractArray(arrayIdx);
  return copyOrCreateCopy(orig, array);
}

}  // namespace

void geoMeshBase::getPointDataArrayCopy(const std::string &arrayName,
                                        vtkAbstractArray *array,
                                        int *arrayIdx) const {
  getArrCopy(_geoMesh.mesh->GetPointData(), arrayName, arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getPointDataArrayCopy(
    const std::string &arrayName, int *arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetPointData(), arrayName, arrayIdx);
}

void geoMeshBase::getPointDataArrayCopy(int arrayIdx,
                                        vtkAbstractArray *array) const {
  getArrCopy(_geoMesh.mesh->GetPointData(), arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getPointDataArrayCopy(
    int arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetPointData(), arrayIdx);
}

void geoMeshBase::getCellDataArrayCopy(const std::string &arrayName,
                                        vtkAbstractArray *array,
                                        int *arrayIdx) const {
  getArrCopy(_geoMesh.mesh->GetCellData(), arrayName, arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getCellDataArrayCopy(
    const std::string &arrayName, int *arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetCellData(), arrayName, arrayIdx);
}

void geoMeshBase::getCellDataArrayCopy(int arrayIdx,
                                        vtkAbstractArray *array) const {
  getArrCopy(_geoMesh.mesh->GetCellData(), arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getCellDataArrayCopy(
    int arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetCellData(), arrayIdx);
}

void geoMeshBase::getFieldDataArrayCopy(const std::string &arrayName,
                                       vtkAbstractArray *array,
                                       int *arrayIdx) const {
  getArrCopy(_geoMesh.mesh->GetFieldData(), arrayName, arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getFieldDataArrayCopy(
    const std::string &arrayName, int *arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetFieldData(), arrayName, arrayIdx);
}

void geoMeshBase::getFieldDataArrayCopy(int arrayIdx,
                                       vtkAbstractArray *array) const {
  getArrCopy(_geoMesh.mesh->GetFieldData(), arrayIdx, array);
}

vtkSmartPointer<vtkAbstractArray> geoMeshBase::getFieldDataArrayCopy(
    int arrayIdx) const {
  return getArrCopy(_geoMesh.mesh->GetFieldData(), arrayIdx);
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
