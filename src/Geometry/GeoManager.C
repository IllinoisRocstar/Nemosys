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
#include "Geometry/GeoManager.H"

#include <unordered_set>

#include <BRepAlgoAPI_BooleanOperation.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRep_Builder.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>

#include "Geometry/ShapeData.H"

namespace NEM {
namespace GEO {

GeoManager::GeoManager(int dim)
    : map_(),
      dim_(std::min(3, std::max(1, dim))) {}

int GeoManager::getDim() const { return dim_; }

void GeoManager::setDim(int dim) { dim_ = dim; }

template <typename Op, typename T, typename F>
std::vector<TopoDS_Shape> modifyTempl(
    Op &&op, T &&shapes, const std::vector<TopAbs_ShapeEnum> &typesToTraverse,
    F &&modifyFunc) {
  std::vector<TopoDS_Shape> shapesToRemove;
  for (auto &shape : shapes) {
    if (typesToTraverse.empty()) {
      modifyFunc(shape, shapesToRemove);
    } else {
      for (auto &shapeType : typesToTraverse) {
        for (TopExp_Explorer explorer{shape, shapeType}; explorer.More();
             explorer.Next()) {
          modifyFunc(explorer.Current(), shapesToRemove);
        }
      }
    }
  }
  return shapesToRemove;
}

std::vector<TopoDS_Shape> GeoManager::modify(
    BRepBuilderAPI_MakeShape &op, const std::vector<TopoDS_Shape> &shapes,
    const std::vector<TopAbs_ShapeEnum> &typesToTraverse) {
  return modifyTempl(op, shapes, typesToTraverse,
                     [this, &op](const TopoDS_Shape &shape,
                                 std::vector<TopoDS_Shape> &shapesToRemove) {
                       this->modifyImpl(op, shape, shapesToRemove);
                     });
}

std::vector<TopoDS_Shape> GeoManager::modify(
    BRepBuilderAPI_MakeShape &op, const TopTools_ListOfShape &shapes,
    const std::vector<TopAbs_ShapeEnum> &typesToTraverse) {
  return modifyTempl(op, shapes, typesToTraverse,
                     [this, &op](const TopoDS_Shape &shape,
                                 std::vector<TopoDS_Shape> &shapesToRemove) {
                       this->modifyImpl(op, shape, shapesToRemove);
                     });
}

std::array<std::vector<TopoDS_Shape>, 2> GeoManager::modify(
    BRepAlgoAPI_BooleanOperation &op,
    const std::vector<TopAbs_ShapeEnum> &typesToTraverse) {
  return {this->modify(op, op.Arguments(), typesToTraverse),
          this->modify(op, op.Tools(), typesToTraverse)};
}

std::vector<TopoDS_Shape> GeoManager::modify(
    BRepBuilderAPI_Sewing &op, const std::vector<TopoDS_Shape> &shapes,
    const std::vector<TopAbs_ShapeEnum> &typesToTraverse) {
  return modifyTempl(
      op, shapes, typesToTraverse,
      [this, &op](const TopoDS_Shape &shape,
                  std::vector<TopoDS_Shape> &shapesToRemove) {
        auto findIter = map_.find(shape);
        if (findIter != map_.end()) {
          auto shapeType = shape.ShapeType();
          auto modified =
              (shapeType == TopAbs_FACE || shapeType == TopAbs_SHELL)
                  ? op.Modified(shape)
                  : op.ModifiedSubShape(shape);
          if (!modified.IsSame(shape)) {
            TopTools_ListOfShape modShapes;
            modShapes.Append(modified);
            if (findIter->second) {
              findIter->second->updateModified(shape, modShapes, *this);
            }
            shapesToRemove.emplace_back(shape);
          }
        }
      });
}

void GeoManager::deleteShapes(const TopoDS_Shape &shape) {
  if (!isChild(shape)) {
    auto findIter = map_.find(shape);
    if (findIter != map_.end()) {
      if (findIter->second) {
        findIter->second->updateDeleted(shape, *this);
      }
      map_.erase(findIter);
    }
  }
}

void GeoManager::deleteShapes(const std::vector<TopoDS_Shape> &shapes) {
  for (auto &shape : shapes) {
    this->deleteShapes(shape);
  }
}

TopoDS_Compound GeoManager::buildCompound() const {
  BRep_Builder builder;
  TopoDS_Compound compound;
  // If no top-dimension shapes, then return a null compound
  bool addedShapes = false;
  for (auto &shape : map_) {
    if (!isChild(shape.first)) {
      if (!addedShapes) {
        builder.MakeCompound(compound);
        addedShapes = true;
      }
      builder.Add(compound, shape.first);
    }
  }
  return compound;
}

std::shared_ptr<ShapeData> *GeoManager::get(const TopoDS_Shape &shape) {
  auto findIter = map_.find(shape);
  if (findIter != map_.end()) {
    return &findIter->second;
  } else {
    return nullptr;
  }
}

const std::shared_ptr<ShapeData> *GeoManager::get(
    const TopoDS_Shape &shape) const {
  auto findIter = map_.find(shape);
  if (findIter != map_.end()) {
    return &findIter->second;
  } else {
    return nullptr;
  }
}

std::pair<GeoManager::MapType::iterator, bool> GeoManager::insert(
    const TopoDS_Shape &shape, std::shared_ptr<ShapeData> shapeData) {
  return map_.emplace(shape, std::move(shapeData));
}

GeoManager::MapType &GeoManager::getMap() { return map_; }

const GeoManager::MapType &GeoManager::getMap() const { return map_; }

bool GeoManager::isChild(const TopoDS_Shape &shape) const {
  int shapeDim;
  switch (shape.ShapeType()) {
    case TopAbs_COMPSOLID:
    case TopAbs_SOLID: shapeDim = 3; break;
    case TopAbs_SHELL:
    case TopAbs_FACE: shapeDim = 2; break;
    case TopAbs_WIRE:
    case TopAbs_EDGE: shapeDim = 1; break;
    case TopAbs_VERTEX: shapeDim = 0; break;
    case TopAbs_COMPOUND:
    case TopAbs_SHAPE:
    default: return false;
  }
  return shapeDim < dim_;
}

void GeoManager::modifyImpl(BRepBuilderAPI_MakeShape &op,
                            const TopoDS_Shape &shape,
                            std::vector<TopoDS_Shape> &shapesToRemove) {
  if (op.IsDeleted(shape)) {
    shapesToRemove.emplace_back(shape);
  } else {
    auto findIter = map_.find(shape);
    if (findIter != map_.end()) {
      auto &generated = op.Generated(shape);
      if (!generated.IsEmpty()) {
        if (findIter->second) {
          findIter->second->updateGenerated(shape, generated, *this);
        }
      }
      auto &modified = op.Modified(shape);
      if (!modified.IsEmpty()) {
        if (findIter->second) {
          findIter->second->updateModified(shape, modified, *this);
        }
        shapesToRemove.emplace_back(shape);
      }
    }
  }
}

}  // namespace NUCMESH
}  // namespace NEM
