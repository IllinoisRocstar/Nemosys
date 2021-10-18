#ifndef NEMOSYS_HEXAGONALARRAY_H_
#define NEMOSYS_HEXAGONALARRAY_H_

#include "nemosys_export.h"
#include "NucMesh/ShapesArray.H"

#include <array>

namespace NEM {
namespace NUCMESH {

class NEMOSYS_EXPORT HexagonalArray : public ShapesArray {
 public:
  HexagonalArray(std::size_t numRadii, double deltaRadius,
                 const std::array<double, 3> &center = {0, 0, 0});

  const std::size_t &getPatternRowCol(int row, int col) const;
  const std::size_t &getPatternCoordCenter(int right, int rightUp) const;

  void setPatternRowCol(int row, int col, std::size_t patternKey);
  /**
   * The shape will be translated by <tt> right * [getGridDistance(), 0, 0] +
   * rightUp * [getGridDistance() * cos(60), getGridDistance() * sin(60), 0] +
   * getCenter() </tt>.
   */
  void setPatternCoordCenter(int right, int rightUp, std::size_t patternKey);

  std::size_t getNumRadii() const { return numRadii_; }
  double getGridDistance() const { return delta_; }
  void setGridDistance(double gridDistance) { delta_ = gridDistance; }

 protected:
  NEM::GEO::GeoManager createGeo() const override;

 private:
  double delta_;
  std::size_t numRadii_;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_HEXAGONALARRAY_H_
