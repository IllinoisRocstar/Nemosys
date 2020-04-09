#ifndef NEMOSYS_SHAPE_H_
#define NEMOSYS_SHAPE_H_

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "nemosys_export.h"

namespace NEM {
namespace GEO {

/**
 * @brief Base class for shapes (circles, polygons, etc.). Provides virtual
 *        methods for shapes.
 */
class NEMOSYS_EXPORT shape {
 public:
  shape() = default;
  virtual ~shape() = default;

 public:
  /**
   * @brief virtual method to draw shapes
   */
  virtual void draw() = 0;

  /**
   * @brief virtual method to apply the mesh type (tri,quad,struct) to surfaces
   */
  virtual void applyMeshType() = 0;

  /**
   * @brief virtual method that maps surfaces to physical entity tag ids
   * @param phystag_map physical entity tags and physical names
   * @param physSurf_map surface and physical tag
   * @return physSurf_map
   */
  virtual std::map<int, int> getPhysSurf(std::map<std::string, int> phystag_map,
                                         std::map<int, int> physSurf_map) = 0;

  /**
   * @brief updates shape surface ids
   * @param
   */
  void updateSurfaces(const std::vector<std::pair<int, int>> &oldNew_vec);

 protected:
  std::vector<double> _center;        /**< center coordinate of shape */
  std::vector<double> _radii;         /**< radii for concentric shapes */
  std::vector<std::string> _meshType; /**< mesh type for each concentric ring */
  std::vector<std::pair<int, int>>
      _elems; /**< number of elements for each concentric ring
               * @note for structured mesh, the pair is (radial,
               * circumferential) directions
               */
  std::vector<std::string> _names; /**< physical entity names (region names) */
  std::vector<int> _surfaces;      /**< container for shape surfaces */
  double _rotation{0.0};           /**< rotation angle (degrees) for polygons */
  int _nSides{3};                  /**< number of side for polygons */
  double _meshSize{0.5};           /**< meshSize at point */
};

}  // namespace GEO
}  // namespace NEM

#endif  // NEMOSYS_SHAPE_H_
