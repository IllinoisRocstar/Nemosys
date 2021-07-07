#ifndef NEMOSYS_CIRCLES_H_
#define NEMOSYS_CIRCLES_H_

#include "nemosys_export.h"
#include "NucMesh/shape.H"

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
   * @param ss bools to make sidesets
   */
  circles(std::vector<double> cen, std::vector<double> rad,
          std::vector<std::string> mtype = {"Tri"},
          std::vector<std::pair<int, int>> el = {{0, 0}},
          std::vector<std::string> name = {},
          std::vector<std::string> ss = {"no"});

  /**
   * @brief Draws or creates the circle
   */
  void draw() override;

  /**
   * @brief Applies the mesh type (tri, quad, struct) to surfaces
   */
  void applyMeshType() override;

  /**
   * @brief Gets all the lines for the circle
   * @return The lines of the circles
   */
  std::vector<std::vector<std::pair<int, int>>> getLines() const;

  /**
   * @brief maps circle surfaces to physical entity tag ids in place
   * @param phystag_map physical entity tags and physical names
   * @param physSurf_map map of surface ids and physical tag ids
   */
  void getPhysSurf(const std::map<std::string, int> &phystag_map,
                   std::map<int, int> &physSurf_map) const override;

  /**
   * @brief Maps circle lines to physical entity tag ids in place
   * @param sstag_map physical entity tags and physical names
   * @param physLine_map map of line ids and physical tag ids
   */
  void getPhysLine(const std::map<std::string, int> &sstag_map,
                   std::map<int, int> &physLine_map) const override;

  /**
   * @brief Returns the max ID/Tag of entity with dimension dim
   * @param dim
   * @return current max id
   */
  int getMaxID(int dim) const;
};

}  // namespace GEO
}  // namespace NEM

#endif  // NEMOSYS_CIRCLES_H_
