#ifndef NEMOSYS_SHAPEBASE_H_
#define NEMOSYS_SHAPEBASE_H_

#include "nemosys_export.h"

#include <array>

namespace NEM {

namespace SRV {
class NucMeshSrv;
}  // namespace SRV

namespace GEO {
class GeoManager;
}  // namespace GEO

namespace NUCMESH {

/**
 * @brief Abstract base class for types that create @c NEM::GEO::GeoManager
 * @details Also see @c NEM::DRV::NucMeshDriver
 */
class NEMOSYS_EXPORT ShapeBase {
 public:
  virtual ~ShapeBase() = default;

  friend class NEM::SRV::NucMeshSrv;
  friend class ShapesArray;

  const std::array<double, 3> &getCenter() const;
  void setCenter(const std::array<double, 3> &center);

  /**
   * @param center Center of circle
   * @param rotation {Radius, angle in degrees}
   * @return Point on circle with center @p center and normal [0, 0, 1], radius
   * @c rotation[0] and angle (in degrees) @c rotation[1]
   */
  static std::array<double, 3> getRotatedPoint(
      const std::array<double, 3> &center,
      const std::array<double, 2> &rotation);

 protected:
  explicit ShapeBase(const std::array<double, 3> &center = {0, 0, 0});
  ShapeBase(const ShapeBase &) = default;
  ShapeBase &operator=(const ShapeBase &) = default;
  ShapeBase(ShapeBase &&) = default;
  ShapeBase &operator=(ShapeBase &&) = default;

  /**
   * Construct a @c NEM::GEO::GeoManager
   * @return A @c NEM::GEO::GeoManager object representing the shape described
   * by this object
   */
  virtual NEM::GEO::GeoManager createGeo() const = 0;

  /**
   * Merge two 2d @c NEM::GEO::GeoManager objects by cutting shapes in @p
   * removeGeo from shapes in @p keepGeo and sewing the resulting faces
   * @param[in,out] keepGeo Merge geometry from @p removeGeo into this object
   * @param[in] removeGeo Move geometry and data from this object into @p
   * keepGeo
   */
  static void mergeGeo(NEM::GEO::GeoManager &keepGeo,
                       NEM::GEO::GeoManager &&removeGeo);

 private:
  std::array<double, 3> center_;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_SHAPEBASE_H_
