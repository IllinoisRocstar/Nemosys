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
   * @param bools to make sidesets
   * @param rotation angle (degrees)
   */
  polygon(int nsides, std::vector<double> cen, std::vector<double> rad,
          std::vector<std::string> mtype = {"Tri"},
          std::vector<std::pair<int, int>> el = {{0, 0}},
          std::vector<std::string> name = {},
          std::vector<std::string> ss = {"no"}, double rot = 0.0);

  /**
   * @brief Draws or creates the polygon
   */
  void draw() override;

  /**
   * @brief Applies the mesh type (tri, quad, struct) to surfaces
   */
  void applyMeshType() override;

  /**
   * @brief Gets all the lines for the polygon
   * @return Vector of all lines of polygon
   */
  std::vector<std::vector<std::pair<int, int>>> getLines() const;

  /**
   * @brief Maps polygon surfaces to physical entity tag ids in place
   * @param phystag_map physical entity tags and physical names
   * @param physSurf_map map of surface ids and physical tag ids
   */
  void getPhysSurf(const std::map<std::string, int> &phystag_map,
                   std::map<int, int> &physSurf_map) const override;

  /**
   * @brief Maps polygon lines to physical entity tag ids in place
   * @param sstag_map physical entity tags and physical names
   * @param physLine_map line ids and physical tag ids
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

#endif  // NEMOSYS_POLYGON_H_
