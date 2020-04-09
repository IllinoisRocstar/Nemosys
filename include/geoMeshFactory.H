#ifndef NEMOSYS_GEOMESHFACTORY_H_
#define NEMOSYS_GEOMESHFACTORY_H_

#include "geoMeshBase.H"

namespace NEM {
namespace MSH {

enum NEMOSYS_EXPORT MeshType { VTK_GEO_MESH, GMSH_GEO_MESH, OSH_GEO_MESH };

NEMOSYS_EXPORT geoMeshBase *Read(const std::string &fileName, MeshType type);
NEMOSYS_EXPORT geoMeshBase *Read(const std::string &fileName);

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_GEOMESHFACTORY_H_
