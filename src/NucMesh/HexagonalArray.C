#include "NucMesh/HexagonalArray.H"

#define _USE_MATH_DEFINES
#include <cmath>

#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

namespace NEM {
namespace NUCMESH {

namespace {

std::size_t coordsToSpiralFlat(int right, int rightUp) {
  /*
   * This converts (right, rightUp) coordinates to a flattened spiral.
   *
   * The value returned is the number of coordinates "dealt with" once
   * that coordinate has been reached, if you start at the center and
   * spiral outwards counterclockwise.
   */
  int ring;     // 0-indexed ring number, counting from the center
  int corner;   // reference corner from which offset is figured
  int offset;   // how far into a side you've gone
  if (right >= 0) {
    if (rightUp >= 0) {
      corner = 0; // start at corner 0 (center-right)
      ring = right + rightUp;
      offset = rightUp; // count upwards
      if (ring == 0) offset = -1; // special case for middle point
    } else if (right <= -rightUp) {
      corner = 4; // start at corner 4 (bottom-left)
      ring = -rightUp;
      offset = right; //count rightwards
    } else {
      corner = 0; // start at corner 0 (center-right)
      ring = right + 1;
      offset = rightUp; // count downwards
    }
  } else {
    if (rightUp <= 0) {
      corner = 3; // start from corner 3 (center-left)
      ring = -(right + rightUp);
      offset = -rightUp; // count (negative) downward
    } else if (-right <= rightUp) {
      corner = 1; // start from corner 1 (top-right)
      ring = rightUp;
      offset = -right; // count (negative) rightwards
    } else {
      corner = 3; // start from corner 3 (center-left)
      ring = -right;
      offset = -rightUp; // count (negative) downward
    }
  }
  // note that 3*N(N-1)+1 is the Nth hexagonal number
  return 3 * ring * (ring - 1) + 1 + corner * ring + offset;
}

void incrementSpiral(int &right, int &rightUp, int &side) {
  switch (side) {
    case 0: {
      --right;
      ++rightUp;
      if (right == 0) ++side;
      break;
    }
    case 1: {
      --right;
      if (right == -rightUp) ++side;
      break;
    }
    case 2: {
      --rightUp;
      if (rightUp == 0) ++side;
      break;
    }
    case 3: {
      ++right;
      --rightUp;
      if (right == 0) ++side;
      break;
    }
    case 4: {
      ++right;
      if (right == -rightUp) ++side;
      break;
    }
    case 5: {
      if (rightUp >= -1) {
        side = 0;
        ++right;
        rightUp = 0;
      } else {
        ++rightUp;
      }
      break;
    }
    default: break;
  }
}

std::array<int, 2> rowColToCoords(int row, int col, std::size_t width) {
  int iwidth = static_cast<int>(width);
  if (row <= iwidth - 1) {
    return {col - iwidth + 1, iwidth - row - 1};
  } else {
    return {row + col - 2 * (iwidth - 1), iwidth - row - 1};
  }
}

}  // namespace

HexagonalArray::HexagonalArray(std::size_t numRadii, double deltaRadius,
                               const std::array<double, 3> &center)
    : ShapesArray(center, 3 * (numRadii) * (numRadii - 1) + 1),
      delta_(deltaRadius),
      numRadii_(numRadii) {}

const std::size_t &HexagonalArray::getPatternRowCol(int row, int col) const {
  auto coords = rowColToCoords(row, col, numRadii_);
  return getPatternCoordCenter(coords[0], coords[1]);
}

const std::size_t &HexagonalArray::getPatternCoordCenter(int right,
                                                         int rightUp) const {
  return this->getPattern(coordsToSpiralFlat(right, rightUp));
}

void HexagonalArray::setPatternRowCol(int row, int col,
                                      std::size_t patternKey) {
  auto coords = rowColToCoords(row, col, numRadii_);
  setPatternCoordCenter(coords[0], coords[1], patternKey);
}

void HexagonalArray::setPatternCoordCenter(int right, int rightUp,
                                           std::size_t patternKey) {
  this->setPattern(coordsToSpiralFlat(right, rightUp), patternKey);
}

// Internal representation of pattern: flattened with all shapes from ring i
// followed by ring i + 1, etc, with shapes in each ring represented
// counterclockwise from positive x-axis
NEM::GEO::GeoManager HexagonalArray::createGeo() const {
  int right = 0;
  int rightUp = 0;
  int side = 5;
  return this->createGeoImpl([this, &right, &rightUp,
                              &side](NEM::GEO::GeoManager *const inp) {
    if (inp) {
      const auto &center = this->getCenter();
      const gp_Vec dest{
          center[0] + delta_ * (right + rightUp * std::cos(60 * M_PI / 180.)),
          center[1] + rightUp * delta_ * std::sin(60 * M_PI / 180.), center[2]};
      gp_Trsf translation{};
      translation.SetTranslation(dest);
      *inp = ShapesArray::basicTransformation(translation, std::move(*inp));
    }
    incrementSpiral(right, rightUp, side);
  });
}

}  // namespace NUCMESH
}  // namespace NEM