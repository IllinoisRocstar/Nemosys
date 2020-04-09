#include "geoMeshFactory.H"

#include "AuxiliaryFunctions.H"

#include "gmshGeoMesh.H"
#ifdef HAVE_OMEGAH
#  include "oshGeoMesh.H"
#endif
#include "vtkGeoMesh.H"

namespace NEM {
namespace MSH {

geoMeshBase *Read(const std::string &fileName) {
  std::string fileExt = nemAux::find_ext(fileName);

  if (fileExt == ".vtu" || fileExt == ".vtk") {
    return Read(fileName, VTK_GEO_MESH);
  } else if (fileExt == ".msh" || fileExt == ".geo") {
    return Read(fileName, GMSH_GEO_MESH);
  } else if (fileExt == ".osh") {
    return Read(fileName, OSH_GEO_MESH);
  } else {
    std::cerr << __func__ << ": File extension " << fileExt
              << " is not supported." << std::endl;
    exit(1);
  }
}

geoMeshBase *Read(const std::string &fileName, MeshType meshType) {
  switch (meshType) {
    case VTK_GEO_MESH: return vtkGeoMesh::Read(fileName);
    case GMSH_GEO_MESH: return gmshGeoMesh::Read(fileName);
    case OSH_GEO_MESH:
#ifdef HAVE_OMEGAH
      return oshGeoMesh::Read(fileName);
#else
      std::cerr << "ERROR: Nemosys is not compiled with Omega_h support."
                << std::endl;
      exit(1);
#endif
  }
  return nullptr;  // should never reach here
}

}  // namespace MSH
}  // namespace NEM
