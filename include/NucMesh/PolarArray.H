#ifndef NEMOSYS_POLARARRAY_H_
#define NEMOSYS_POLARARRAY_H_

#include "nemosys_export.h"
#include "NucMesh/ShapesArray.H"

#include <array>

namespace NEM {
namespace NUCMESH {

/**
 * Arrange a series of other @c ShapeBase objects along a circular arc.
 */
class NEMOSYS_EXPORT PolarArray : public ShapesArray {
 public:
  PolarArray(std::size_t numSubshapes, double startAngle, double endAngle,
             double radius, bool rotateWithArray = false,
             const std::array<double, 3> &center = {0, 0, 0});

  using ShapesArray::getPattern;
  /**
   * The shape at pattern index @c i will be, for <tt> theta = i *
   * std::fmod(getEndAngle() - getStartAngle(), 360.) </tt>, translated by <tt>
   * ShapeBase::getRotatedPoint(getCenter(), {getRadius(), theta}) </tt> and (if
   * @c rotatewithArray) rotated about [0, 0, 1] by @c theta.
   */
  using ShapesArray::setPattern;

  std::size_t getNumSubshapes() const { return numShapesInArr_; }
  double getRadius() const { return radius_; }
  void setRadius(double radius) { radius_ = radius; }
  double getStartAngle() const { return start_; }
  void setStartAngle(double startAngle) { start_ = startAngle; }
  double getEndAngle() const { return end_; }
  void setEndAngle(double endAngle) { end_ = endAngle; }
  bool getRotateWithArray() const { return rotateWithArray_; }
  void setRotateWithArray(bool rotateWithArray) {
    rotateWithArray_ = rotateWithArray;
  }

 protected:
  NEM::GEO::GeoManager createGeo() const override;

 private:
  std::size_t numShapesInArr_;
  double radius_;
  double start_;
  double end_;
  bool rotateWithArray_;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_POLARARRAY_H_
