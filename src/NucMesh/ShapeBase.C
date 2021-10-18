#include "NucMesh/ShapeBase.H"

#define _USE_MATH_DEFINES
#include <cmath>

#include <array>
#include <vector>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_ListOfShape.hxx>
#include <TopoDS_Shape.hxx>

#include "Geometry/GeoManager.H"

namespace NEM {
namespace NUCMESH {

const std::array<double, 3> &ShapeBase::getCenter() const { return center_; }

void ShapeBase::setCenter(const std::array<double, 3> &center) {
  center_ = center;
}

ShapeBase::ShapeBase(const std::array<double, 3> &center) : center_(center) {}

std::array<double, 3> ShapeBase::getRotatedPoint(
    const std::array<double, 3> &center,
    const std::array<double, 2> &rotation) {
  auto angle = rotation[1] * M_PI / 180.;
  return {center[0] + rotation[0] * std::cos(angle),
          center[1] + rotation[0] * std::sin(angle), center[2]};
}

void ShapeBase::mergeGeo(NEM::GEO::GeoManager &keepGeo,
                         NEM::GEO::GeoManager &&removeGeo) {
  auto origCompound = keepGeo.buildCompound();
  if (origCompound.IsNull()) {
    keepGeo = std::move(removeGeo);
    return;
  } else {
    auto newCompound = removeGeo.buildCompound();
    TopoDS_ListOfShape newShapes;
    static constexpr std::array<TopAbs_ShapeEnum, 4> shapeTypes{
        TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX};
    for (auto &shapeType : shapeTypes) {
      for (TopExp_Explorer explorer{newCompound, shapeType}; explorer.More();
           explorer.Next()) {
        auto &oldSubshape = explorer.Current();
        if (auto old_metadata = removeGeo.get(oldSubshape)) {
          keepGeo.insert(oldSubshape, std::move(*old_metadata));
        }
      }
    }
    if (!newCompound.IsNull()) {
      keepGeo.setDim(std::max(keepGeo.getDim(), removeGeo.getDim()));
      BRepAlgoAPI_Cut cutter{origCompound, newCompound};
      if (cutter.HasDeleted() || cutter.HasModified() ||
          cutter.HasGenerated()) {
        auto deletedShapes = keepGeo.modify(cutter, {TopAbs_EDGE, TopAbs_FACE});
        keepGeo.deleteShapes(deletedShapes[0]);
      }
      BRepBuilderAPI_Sewing sewer{};
      std::vector<TopoDS_Shape> sewedShapes;
      auto compound = keepGeo.buildCompound();
      for (TopExp_Explorer explorer{compound, TopAbs_FACE}; explorer.More();
           explorer.Next()) {
        sewer.Add(explorer.Current());
        sewedShapes.emplace_back(explorer.Current());
      }
      sewer.Perform();
      auto deletedShapes =
          keepGeo.modify(sewer, sewedShapes, {TopAbs_EDGE, TopAbs_FACE});
      keepGeo.deleteShapes(deletedShapes);
    }
  }
}

}  // namespace NUCMESH
}  // namespace NEM
