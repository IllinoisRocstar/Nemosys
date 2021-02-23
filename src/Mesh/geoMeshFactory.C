#include "geoMeshFactory.H"

#include "AuxiliaryFunctions.H"

#include "inpGeoMesh.H"
#include "exoGeoMesh.H"
#include "gmshGeoMesh.H"
#include "oshGeoMesh.H"
#include "vtkGeoMesh.H"

namespace NEM {
namespace MSH {

MeshType MeshTypeFromFilename(const std::string &fileName) {
  std::string fileExt = nemAux::find_ext(fileName);

  if (fileExt == ".vtu" || fileExt == ".vtk") {
    return MeshType::VTK_GEO_MESH;
  } else if (fileExt == ".msh" || fileExt == ".geo") {
    return MeshType::GMSH_GEO_MESH;
  } else if (fileExt == ".osh") {
    return MeshType::OSH_GEO_MESH;
  } else if (fileExt == ".exo" || fileExt == ".e" || fileExt == ".gen" ||
             fileExt == ".g") {
    return MeshType::EXO_GEO_MESH;
  } else if (fileExt == ".inp") {
    return MeshType::INP_GEO_MESH;
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
    case MeshType::VTK_GEO_MESH: return vtkGeoMesh::Read(fileName);
    case MeshType::GMSH_GEO_MESH: return gmshGeoMesh::Read(fileName);
    case MeshType::OSH_GEO_MESH: return oshGeoMesh::Read(fileName);
    case MeshType::EXO_GEO_MESH: return exoGeoMesh::Read(fileName);
    case MeshType::INP_GEO_MESH: return inpGeoMesh::Read(fileName);
  }
}

geoMeshBase *New(MeshType meshType) {
  switch (meshType) {
    case MeshType::VTK_GEO_MESH: return vtkGeoMesh::New();
    case MeshType::GMSH_GEO_MESH: return gmshGeoMesh::New();
    case MeshType::OSH_GEO_MESH: return oshGeoMesh::New();
    case MeshType::EXO_GEO_MESH: return exoGeoMesh::New();
    case MeshType::INP_GEO_MESH: return inpGeoMesh::New();
  }
}

geoMeshBase *New(const std::string &fileName) {
  return New(MeshTypeFromFilename(fileName));
}

}  // namespace MSH
}  // namespace NEM
