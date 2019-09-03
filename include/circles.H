#ifndef NEMOSYS_CIRCLES_H_
#define NEMOSYS_CIRCLES_H_

#include "nemosys_export.h"
#include "shape.H"

namespace NEM {
namespace GEO {

/**
 * @brief Class for constructing concentric circles
 */
class NEMOSYS_EXPORT circles : public shape {
 public:
  /**
   * @brief Constructs either a single circle or concentric circles.
   * @param cen center point
   * @param rad radii for each circle
   * @param mtype mesh type for each circle
   * @param el number of elements if structured (radial, circum)
   * @param name region name
   */
  circles(std::vector<double> cen, std::vector<double> rad,
          std::vector<std::string> mtype = {"Tri"},
          std::vector<std::pair<int, int>> el = {{0, 0}},
          std::vector<std::string> name = {});

  /**
   * @brief Draws or creates the circle
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

#endif  // NEMOSYS_CIRCLES_H_
