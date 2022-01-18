/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_SHAPESARRAY_H_
#define NEMOSYS_SHAPESARRAY_H_

#include "nemosys_export.h"
#include "NucMesh/ShapeBase.H"

#include <map>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "Geometry/GeoManager.H"

class gp_Trsf;

namespace NEM {
namespace NUCMESH {

/**
 * Abstract base class representing a set of other @c ShapeBase objects, with a
 * transformation applied to each.
 */
class NEMOSYS_EXPORT ShapesArray : public ShapeBase {
 public:
  std::size_t getNumPatternShapes() const;

  const ShapeBase *getPatternShape(std::size_t idx) const;

  void setPatternShape(std::size_t idx,
                       const std::shared_ptr<ShapeBase> &shape);

  template <typename Shape>
  void insertPatternShape(std::size_t idx, Shape &&shape) {
    this->setPatternShape(idx,
                          std::make_shared<typename std::decay<Shape>::type>(
                              std::forward<Shape>(shape)));
  }

  template <typename Shape, typename... Args>
  void makePatternShape(std::size_t idx, Args &&...args) {
    this->setPatternShape(idx,
                          std::make_shared<Shape>(std::forward<Args>(args)...));
  }

  void fillPattern(std::size_t idx);

 protected:
  ShapesArray(const std::array<double, 3> &center,
              std::size_t numPatternShapes);

  std::size_t getPatternSize() const;

  const std::size_t &getPattern(std::size_t idx) const;

  void setPattern(std::size_t idx, std::size_t patternKey);

  template <typename Modifier>
  NEM::GEO::GeoManager createGeoImpl(Modifier &&modifier) const {
    NEM::GEO::GeoManager output(0);
    for (auto &patternId : pattern_) {
      auto shape = patternId < patternShapes_.size()
                       ? patternShapes_.at(patternId).get()
                       : nullptr;
      if (shape) {
        auto subShapeOut = shape->createGeo();
        modifier(&subShapeOut);
        ShapeBase::mergeGeo(output, std::move(subShapeOut));
      } else {
        modifier(nullptr);
      }
    }
    return output;
  }

  static NEM::GEO::GeoManager basicTransformation(
      const gp_Trsf &transformation, NEM::GEO::GeoManager &&geoMetadata);

 private:
  /**
   * A set of other @c ShapeBase objects, referenced by @c pattern_
   */
  std::vector<std::shared_ptr<ShapeBase>> patternShapes_;
  /**
   * Each entry of @c pattern_ is an index into @c patternShapes_; subclasses
   * interpret how to transform each entry
   */
  std::vector<decltype(patternShapes_)::size_type> pattern_;
};

}  // namespace NUCMESH
}  // namespace NEM

#endif  // NEMOSYS_SHAPESARRAY_H_
