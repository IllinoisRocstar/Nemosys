#include "geoMeshFactory.H"

#include "AuxiliaryFunctions.H"

#include "exoGeoMesh.H"
#include "gmshGeoMesh.H"
#include "oshGeoMesh.H"
#include "vtkGeoMesh.H"

namespace NEM {
namespace MSH {

MeshType MeshTypeFromFilename(const std::string &fileName) {
  std::string fileExt = nemAux::find_ext(fileName);

  if (fileExt == ".vtu" || fileExt == ".vtk") {
    return VTK_GEO_MESH;
  } else if (fileExt == ".msh" || fileExt == ".geo") {
    return GMSH_GEO_MESH;
  } else if (fileExt == ".osh") {
    return OSH_GEO_MESH;
  } else if (fileExt == ".exo" || fileExt == ".e" || fileExt == ".gen" ||
             fileExt == ".g") {
    return EXO_GEO_MESH;
  } else {
    std::cerr << "File extension " << fileExt << " is not supported."
              << std::endl;
    exit(1);
  }
}

geoMeshBase *Read(const std::string &fileName) {
  return Read(fileName, MeshTypeFromFilename(fileName));
}

geoMeshBase *Read(const std::string &fileName, MeshType meshType) {
  switch (meshType) {
    case VTK_GEO_MESH: return vtkGeoMesh::Read(fileName);
    case GMSH_GEO_MESH: return gmshGeoMesh::Read(fileName);
    case OSH_GEO_MESH: return oshGeoMesh::Read(fileName);
    case EXO_GEO_MESH: return exoGeoMesh::Read(fileName);
  }
  return nullptr;  // should never reach here
}

geoMeshBase *New(MeshType meshType) {
  switch (meshType) {
    case VTK_GEO_MESH: return vtkGeoMesh::New();
    case GMSH_GEO_MESH: return gmshGeoMesh::New();
    case OSH_GEO_MESH: return oshGeoMesh::New();
    case EXO_GEO_MESH: return exoGeoMesh::New();
  }
  return nullptr;  // should never reach here
}

geoMeshBase *New(const std::string &fileName) {
  return New(MeshTypeFromFilename(fileName));
}

}  // namespace MSH
}  // namespace NEM
