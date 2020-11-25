#ifndef NEMOSYS_POLYGON_H_
#define NEMOSYS_POLYGON_H_

#include "nemosys_export.h"
#include "shape.H"

namespace NEM {
namespace GEO {

/**
 * @brief Class for constructing concentric polygons
 */
class NEMOSYS_EXPORT polygon : public shape {
 public:
  /**
   * @brief Constructs either a single polygon or concentric polygons.
   * @param center point
   * @param radii of points
   * @param mesh type
   * @param number of elements if structured (radial, circum)
   * @param region name
   * @param rotation angle (degrees)
   */
  polygon(int nsides, std::vector<double> cen, std::vector<double> rad,
          std::vector<std::string> mtype = {"Tri"},
          std::vector<std::pair<int, int>> el = {{0, 0}},
          std::vector<std::string> name = {}, double rot = 0.0);

  /**
   * @brief Draws or creates the polygon
   */
  void draw() override;

  /**
   * @brief Applies the mesh type (tri, quad, struct) to surfaces
   */
  void applyMeshType() override;

  /**
   * @brief maps circle surfaces to physical entity tag ids
   * @param phystag_map physical entity tags and physical names
   * @param physSurf_map surface and physical tag
   * @return physSurf_map
   */
  std::map<int, int> getPhysSurf(std::map<std::string, int> phystag_map,
                                 std::map<int, int> physSurf_map) override;

  /**
   * @brief Returns the max ID/Tag of entity with dimension dim
   * @param dim
   * @return current max id
   */
  static int getMaxID(int dim);
};

}  // namespace GEO
}  // namespace NEM

#endif  // NEMOSYS_POLYGON_H_
