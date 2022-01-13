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
#ifndef NEMOSYS_GEOMANAGER_H_
#define NEMOSYS_GEOMANAGER_H_

#include "nemosys_export.h"

#include <array>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <TopAbs_ShapeEnum.hxx>
#include <TopTools_ShapeMapHasher.hxx>

#include "Geometry/ShapeData.H"

class TopoDS_Compound;
class TopoDS_Shape;
class BRepBuilderAPI_MakeShape;
class BRepAlgoAPI_BooleanOperation;
class BRepBuilderAPI_Sewing;

namespace NEM {
namespace GEO {

/**
 * @brief Class to manage @c TopoDS_Shapes along with metadata
 * @details Handles results of @c BRepBuilderAPI_MakeShape (in particular,
 * boolean operations) by delegating to the arguments' metadata. Also see @c
 * NEM::GEO::ShapeData
 */
class NEMOSYS_EXPORT GeoManager {
 public:
  /**
   * Hasher based on TopTools_ShapeMapHasher, which relies on
   * TopoDS_Shape::IsSame, which is not TopoDS_Shape::operator==
   */
  struct ShapeMapHasher_Hash {
    std::size_t operator()(const TopoDS_Shape &shape) const {
      return TopTools_ShapeMapHasher::HashCode(shape, IntegerLast());
    }
  };

  /**
   * Equality operator based on TopTools_ShapeMapHasher, which relies on
   * TopoDS_Shape::IsSame, which is not TopoDS_Shape::operator==
   */
  struct ShapeMapHasher_KeyEqual {
    bool operator()(const TopoDS_Shape &shape1,
                    const TopoDS_Shape &shape2) const {
      return TopTools_ShapeMapHasher::IsEqual(shape1, shape2);
    }
  };
  using MapType =
      std::unordered_map<TopoDS_Shape, std::shared_ptr<ShapeData>,
                         ShapeMapHasher_Hash, ShapeMapHasher_KeyEqual>;

  /**
   * Create an empty geometry manager
   * @param dim Shapes of this dimension (clipped to be between 1 and 3,
   * inclusive) are treated as first-class citizens, with lower-dimensional
   * shapes assumed to be subshapes.
   */
  explicit GeoManager(int dim);

  /**
   * Get dimension of geometry
   * @return dimension
   */
  int getDim() const;
  /**
   * Set the dimension
   * @param dim dimension (between 1 and 3)
   */
  void setDim(int dim);

  /**
   * Modify this geoMetadata after a @c BRepBuilderAPI_MakeShape operation
   * @param op Modification of geometry
   * @param shapes Geometry objects assumed to be inputs of @p op. The
   * corresponding values in the map dictate the behavior of the resulting
   * metadata.
   * @param typesToTraverse Traverse the subshapes of @p shapes with these
   * types. An empty list means only @p shapes will be examined.
   * @return Shapes that are no longer present in the output of @p op
   */
  std::vector<TopoDS_Shape> modify(
      BRepBuilderAPI_MakeShape &op, const std::vector<TopoDS_Shape> &shapes,
      const std::vector<TopAbs_ShapeEnum> &typesToTraverse = {
          TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID});

  /**
   * @copydoc GeoManager::modify(BRepBuilderAPI_MakeShape &, const
   * std::vector<TopoDS_Shape> &, const std::vector<TopAbs_ShapeEnum> &)
   */
  std::vector<TopoDS_Shape> modify(
      BRepBuilderAPI_MakeShape &op, const TopTools_ListOfShape &shapes,
      const std::vector<TopAbs_ShapeEnum> &typesToTraverse = {
          TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID});

  /**
   * @copydoc GeoManager::modify(BRepBuilderAPI_MakeShape &, const
   * std::vector<TopoDS_Shape> &, const std::vector<TopAbs_ShapeEnum> &)
   * shapes are inferred from @op
   */
  std::array<std::vector<TopoDS_Shape>, 2> modify(
      BRepAlgoAPI_BooleanOperation &op,
      const std::vector<TopAbs_ShapeEnum> &typesToTraverse = {
          TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID});

  /**
   * @copydoc GeoManager::modify(BRepBuilderAPI_MakeShape &, const
   * std::vector<TopoDS_Shape> &, const std::vector<TopAbs_ShapeEnum> &)
   */
  std::vector<TopoDS_Shape> modify(
      BRepBuilderAPI_Sewing &op, const std::vector<TopoDS_Shape> &shapes,
      const std::vector<TopAbs_ShapeEnum> &typesToTraverse = {
          TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID});

  /**
   * @brief Remove a shape from the map
   * @details Note there is no effect if shape's dimension is less than the
   * instance dimension (to avoid removing shapes that are still referenced by
   * other shapes).
   * @param shape Shape to remove, if present
   */
  void deleteShapes(const TopoDS_Shape &shape);

  /**
   * @brief Delete a set of shapes
   * @copydetails GeoManager::deleteShapes(const TopoDS_Shape &)
   * @param shapes Shapes to remove, if present
   */
  void deleteShapes(const std::vector<TopoDS_Shape> &shapes);

  /**
   * Create a compound from shapes present in the map that have same dimension
   * as the instance.
   * @return compound of shapes
   */
  TopoDS_Compound buildCompound() const;

  /**
   * @copydoc GeoManager::get(const TopoDS_Shape &) const
   */
  std::shared_ptr<ShapeData> *get(const TopoDS_Shape &shape);

  /**
   * Get the data from the map
   * @param shape Shape to query
   * @return Pointer to data in the map, or nullptr if shape not in map. If a
   * pointer to a nullptr, then shape is in map, but data in map is null.
   */
  const std::shared_ptr<ShapeData> *get(const TopoDS_Shape &shape) const;

  std::pair<MapType::iterator, bool> insert(
      const TopoDS_Shape &shape, std::shared_ptr<ShapeData> shapeData);

  /**
   * Helper to insert a shape and data into the map
   * @tparam T @c ShapeData type
   * @param shape Shape to insert
   * @param shapeData Data to associate to with @shape
   * @return iterator to inserted or already existing element, and whether or
   * not there was an insertion
   */
  template <typename T>
  std::pair<MapType::iterator, bool> insertForward(const TopoDS_Shape &shape,
                                                   T &&shapeData) {
    auto emplaceIter = getMap().emplace(shape, std::shared_ptr<ShapeData>{});
    if (emplaceIter.second) {
      emplaceIter.first->second.reset(
          new typename std::decay<T>::type(std::forward<T>(shapeData)));
    }
    return emplaceIter;
  }

  /**
   * Helper to insert a shape and construct data
   * @tparam T @c ShapeData type
   * @tparam Args Types of @p T constructor
   * @param shape Shape to insert
   * @param args Arguments of @p T constructor for data to associate with @p
   * shape
   * @return iterator to inserted or already existing element, and whether or
   * not there was an insertion
   */
  template <typename T, typename... Args>
  std::pair<MapType::iterator, bool> insertConstruct(const TopoDS_Shape &shape,
                                                     Args &&...args) {
    auto emplaceIter = getMap().emplace(shape, std::shared_ptr<ShapeData>{});
    if (emplaceIter.second) {
      emplaceIter.first->second.reset(new T(std::forward<Args>(args)...));
    }
    return emplaceIter;
  }

  /*
   * @copydoc GeoManager::getMap() const
   */
  MapType &getMap();

  /**
   * Get the shapes and data
   * @return Map from shapes to data
   */
  const MapType &getMap() const;

 protected:
  MapType map_;
  int dim_;

  bool isChild(const TopoDS_Shape &shape) const;

  void modifyImpl(BRepBuilderAPI_MakeShape &op, const TopoDS_Shape &shape,
                  std::vector<TopoDS_Shape> &shapesToRemove);
};

}  // namespace GEO
}  // namespace NEM

#endif  // NEMOSYS_GEOMANAGER_H_
