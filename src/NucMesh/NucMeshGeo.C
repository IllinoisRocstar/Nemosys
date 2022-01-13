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
#include "NucMesh/NucMeshGeo.H"

#include <unordered_set>

#include <TopExp_Explorer.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Shape.hxx>

#include <SMESH_Gen.hxx>
#include <SMESH_Mesh.hxx>
#include <StdMeshers_Adaptive1D.hxx>
#include <StdMeshers_MEFISTO_2D.hxx>
#include <StdMeshers_Regular_1D.hxx>

#include "NucMesh/NucMeshShapeData.H"

namespace NEM {
namespace NUCMESH {

std::unique_ptr<SMESH_Mesh> NucMeshGeo::computeMesh(SMESH_Gen &generator) {
  std::unique_ptr<SMESH_Mesh> mesh{generator.CreateMesh(false)};
  std::vector<std::unique_ptr<SMESH_Hypothesis>> hypotheses;
  auto compound = this->buildCompound();
  mesh->ShapeToMesh(compound);
  {
    auto edgeAlgId = generator.GetANewId();
    hypotheses.emplace_back(new StdMeshers_Regular_1D{edgeAlgId, &generator});
    mesh->AddHypothesis(compound, edgeAlgId);
    auto edgeHypId = generator.GetANewId();
    hypotheses.emplace_back(
        new StdMeshers_Adaptive1D{edgeHypId, &generator});
    mesh->AddHypothesis(compound, edgeHypId);
  }
  {
    auto faceAlgId = generator.GetANewId();
    hypotheses.emplace_back(new StdMeshers_MEFISTO_2D{faceAlgId, &generator});
    mesh->AddHypothesis(compound, faceAlgId);
  }
  std::unordered_set<TopoDS_Shape, ShapeMapHasher_Hash, ShapeMapHasher_KeyEqual>
      generated_shapes;
  for (auto &shape : map_) {
    if (!isChild(shape.first)) {
      static constexpr std::array<TopAbs_ShapeEnum, 4> shapeTypes{
          TopAbs_VERTEX, TopAbs_EDGE, TopAbs_FACE, TopAbs_SOLID};
      for (auto &shapeType : shapeTypes) {
        for (TopExp_Explorer explorer{shape.first, shapeType}; explorer.More();
             explorer.Next()) {
          auto &subShape = explorer.Current();
          auto findIter = map_.find(subShape);
          if (findIter != map_.end() &&
              generated_shapes.find(subShape) == generated_shapes.end()) {
            if (auto nmData =
                    dynamic_cast<NucMeshShapeData *>(findIter->second.get())) {
              nmData->setupAlgos(explorer.Current(), generator, *mesh,
                                 hypotheses);
            }
            generated_shapes.emplace(subShape);
          }
        }
      }
    }
  }
  generator.Compute(*mesh, compound);
  return mesh;
}

}  // namespace NUCMESH
}  // namespace NEM