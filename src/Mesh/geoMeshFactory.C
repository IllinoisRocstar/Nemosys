#include "Mesh/geoMeshFactory.H"

#include "AuxiliaryFunctions.H"

#include "Mesh/exoGeoMesh.H"
#include "Mesh/inpGeoMesh.H"
#include "Mesh/oshGeoMesh.H"
#include "Mesh/vtkGeoMesh.H"

#ifdef HAVE_CFMSH
#  include "Mesh/foamGeoMesh.H"
#endif

#ifdef HAVE_GMSH
#  include "Mesh/gmshGeoMesh.H"
#endif

#ifdef HAVE_OCC
#  include "Mesh/smeshGeoMesh.H"
#endif

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
  } else if (fileExt == ".foam") {
    return MeshType::FOAM_GEO_MESH;
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
    case MeshType::GMSH_GEO_MESH: {
#ifdef HAVE_GMSH
      return gmshGeoMesh::Read(fileName);
#else
      std::cerr << "Please build NEMoSys with ENABLE_GMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::OSH_GEO_MESH: return oshGeoMesh::Read(fileName);
    case MeshType::EXO_GEO_MESH: return exoGeoMesh::Read(fileName);
    case MeshType::INP_GEO_MESH: return inpGeoMesh::Read(fileName);
    case MeshType::FOAM_GEO_MESH: {
#ifdef HAVE_CFMSH
      std::string fname = fileName;
      fname.erase(fname.find_last_of('.'));
      return foamGeoMesh::Read(fname);
#else
      std::cerr << "Please build NEMoSys with ENABLE_CFMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::SMESH_GEO_MESH: {
      std::cerr << "Unsupported.\n";
      throw;
    }
  }
  // Don't put inside switch statement so static analyzers can detect unhandled
  // cases
  return nullptr;
}

geoMeshBase *New(MeshType meshType) {
  switch (meshType) {
    case MeshType::VTK_GEO_MESH: return vtkGeoMesh::New();
    case MeshType::GMSH_GEO_MESH: {
#ifdef HAVE_GMSH
      return gmshGeoMesh::New();
#else
      std::cerr << "Please build NEMoSys with ENABLE_GMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::OSH_GEO_MESH: return oshGeoMesh::New();
    case MeshType::EXO_GEO_MESH: return exoGeoMesh::New();
    case MeshType::INP_GEO_MESH: return inpGeoMesh::New();
    case MeshType::FOAM_GEO_MESH: {
#ifdef HAVE_CFMSH
      return foamGeoMesh::New();
#else
      std::cerr << "Please build NEMoSys with ENABLE_CFMSH=ON to use this"
                << " feature!" << std::endl;
      throw;
#endif
    }
    case MeshType::SMESH_GEO_MESH: {
#ifdef HAVE_OCC
      return smeshGeoMesh::New();
#else
      std::cerr << "Requires ENABLE_OPENCASCADE=ON\n";
      throw;
#endif
    }
  }
  return nullptr;
}

geoMeshBase *New(const std::string &fileName) {
  return New(MeshTypeFromFilename(fileName));
}

}  // namespace MSH
}  // namespace NEM
