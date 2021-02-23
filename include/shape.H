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
   * @brief virtual method that maps surfaces, in place, to physical entity tag
   * ids
   * @param phystag_map physical entity tags and physical names
   * @param physSurf_map map containing surface id and physical tag id
   */
  virtual void getPhysSurf(const std::map<std::string, int> &phystag_map,
                           std::map<int, int> &physSurf_map) const = 0;

  /**
   * @brief virtual method that maps lines to physical entity tag ids in place
   * @param sstag_map physical entity tags and physical names
   * @param physLine_map line ids and physical tag ids
   */
  virtual void getPhysLine(const std::map<std::string, int> &sstag_map,
                           std::map<int, int> &physLine_map) const = 0;

  /**
   * @brief updates shape surface ids
   * @param oldNew_vec vector of pairs of old-to-new surface ids
   */
  void updateSurfaces(const std::vector<std::pair<int, int>> &oldNew_vec);

 protected:
  std::vector<double> _center;        /**< center coordinate of shape */
  std::vector<double> _radii;         /**< radii for concentric shapes */
  std::vector<double> _circle_radii;  /**< radii for concentric circle */
  std::vector<double> _poly_radii;    /**< radii for concentric polygon */
  std::vector<std::string> _meshType; /**< mesh type for each concentric ring */
  std::vector<std::pair<int, int>>
      _elems; /**< number of elements for each concentric ring
               * @note for structured mesh, the pair is (radial,
               * circumferential) directions
               */
  std::vector<std::string> _names; /**< physical entity names (region names) */
  std::vector<std::string> _sideset; /**< bool to make sideset */
  std::vector<int> _surfaces;        /**< container for shape surfaces */
  double _rotation{0.0}; /**< rotation angle (degrees) for polygons */
  int _nSides{3};        /**< number of sides for polygons */
  double _meshSize{0.5}; /**< meshSize at point */
};

}  // namespace GEO
}  // namespace NEM

#endif  // NEMOSYS_SHAPE_H_
