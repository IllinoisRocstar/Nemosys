#ifndef NEMOSYS_CIRCLESINPOLYS_H_
#define NEMOSYS_CIRCLESINPOLYS_H_

#include "nemosys_export.h"
#include "shape.H"

namespace NEM {
namespace GEO {

/**
 * @brief Class for constructing concentric polygons
 */
class NEMOSYS_EXPORT circlesInPolys : public shape {
 public:
  /**
   * @brief Constructs concentric circles inside of a polygon.
   * @param center point
   * @param radii of circle points
   * @param radii of polygon points
   * @param mesh type
   * @param number of elements if structured (radial, circum)
   * @param region name
   * @param sidesets bools
   * @param rotation angle (degrees)
   */
  circlesInPolys(int nsides, std::vector<double> cen,
                 std::vector<double> circle_rad, std::vector<double> poly_rad,
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
   * @brief Gets all the lines for the circles in polygons
   * @return Vector of the lines of the circle in polygon
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
   * @brief maps circlesInPolys lines to physical entity tag ids in place
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

#endif  // NEMOSYS_CIRCLESINPOLYS_H_
