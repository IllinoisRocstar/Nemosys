#include "NucMesh/PolarArray.H"

#define _USE_MATH_DEFINES
#include <cmath>

#include <gp_Ax1.hxx>
#include <gp_Trsf.hxx>

namespace NEM {
namespace NUCMESH {

namespace {

double fixEndAngle(double startAngle, double endAngle) {
  auto mapped = std::fmod(endAngle, 360.);
  return mapped > startAngle ? mapped : mapped + 360;
}

}  // namespace

PolarArray::PolarArray(std::size_t numSubshapes, double startAngle,
                       double endAngle, double radius, bool rotateWithArray,
                       const std::array<double, 3> &center)
    : ShapesArray(center, numSubshapes),
      numShapesInArr_(numSubshapes),
      radius_(radius),
      start_(startAngle),
      end_(endAngle),
      rotateWithArray_(rotateWithArray) {}

NEM::GEO::GeoManager PolarArray::createGeo() const {
  const auto start = std::fmod(start_, 360.);
  const auto end = fixEndAngle(start, end_);
  const auto increment = (end - start) / numShapesInArr_;
  double angle = start;
  return this->createGeoImpl([this, increment,
                              &angle](NEM::GEO::GeoManager *const inp) {
    if (inp) {
      const auto destCenter =
          getRotatedPoint(this->getCenter(), {this->radius_, angle});
      gp_Trsf transformation{};
      transformation.SetTranslation(
          gp_Vec{destCenter[0], destCenter[1], destCenter[2]});
      if (this->rotateWithArray_) {
        gp_Trsf rotation{};
        gp_Ax1 axisOfRotation{};  // Z-axis by default
        rotation.SetRotation(axisOfRotation, angle * M_PI / 180.);
        transformation *= rotation;
      }
      *inp = ShapesArray::basicTransformation(transformation, std::move(*inp));
    }
    angle += increment;
  });
}

}  // namespace NUCMESH
}  // namespace NEM