#ifndef NEMOSYS_NUCMESHGEO_H_
#define NEMOSYS_NUCMESHGEO_H_

#include "nemosys_export.h"
#include "Geometry/GeoManager.H"

#include <memory>
#include <utility>

class SMESH_Mesh;
class SMESH_Gen;

namespace NEM {
namespace NUCMESH {

class NEMOSYS_EXPORT NucMeshGeo : public NEM::GEO::GeoManager {
  using NEM::GEO::GeoManager::GeoManager;
 public:
  /**
   * @brief Construct a mesh using the shapes in the map and their @c ShapeData.
   * @details Note there is a default edge hypothesis
   * (StdMeshers_AutomaticLength), edge algorithm (StdMeshers_Regular_1D), and
   * face algorithm (StdMeshers_MEFISTO_2D). Geometry and mesh regenerated on
   * each call.
   * @return Generated mesh
   */
  std::unique_ptr<SMESH_Mesh> computeMesh(SMESH_Gen &generator);
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_NUCMESHGEO_H_
