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

#include <NETGENPlugin_Hypothesis_2D.hxx>
#include <NETGENPlugin_NETGEN_2D.hxx>
#include <NETGENPlugin_NETGEN_2D_ONLY.hxx>
#include <SMESH_Gen.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESH_MeshEditor.hxx>
#include <StdMeshers_Adaptive1D.hxx>
#include <StdMeshers_MEFISTO_2D.hxx>
#include <StdMeshers_Regular_1D.hxx>

#include "NucMesh/NucMeshShapeData.H"

namespace NEM {
namespace NUCMESH {

std::unique_ptr<SMESH_Mesh> NucMeshGeo::computeMesh(SMESH_Gen &generator,
                                                    bool useNetgen,
                                                    double maxMeshSize,
                                                    double minMeshSize) {
  algosTime_.start();
  std::unique_ptr<SMESH_Mesh> mesh{generator.CreateMesh(false)};
  std::vector<std::unique_ptr<SMESH_Hypothesis>> hypotheses;
  auto compound = this->buildCompound();
  mesh->ShapeToMesh(compound);
  {
    auto edgeAlgId = generator.GetANewId();
    // Set Edge algo
    hypotheses.emplace_back(new StdMeshers_Regular_1D{edgeAlgId, &generator});
    mesh->AddHypothesis(compound, edgeAlgId);
    auto edgeHypId = generator.GetANewId();
    // Set Edge mesh hypothesis
    auto edgeHypo = new StdMeshers_Adaptive1D{edgeHypId, &generator};
    // Set the mesh sizes
    if (minMeshSize < 1.e12) {
      edgeHypo->SetMinSize(minMeshSize);
      edgeHypo->SetMaxSize(maxMeshSize);
      edgeHypo->SetDeflection(minMeshSize);
    }
    hypotheses.emplace_back(edgeHypo);
    mesh->AddHypothesis(compound, edgeHypId);
  }
  {
    auto faceAlgId = generator.GetANewId();
    if (useNetgen) {
#ifdef HAVE_NGEN
      auto alg = new NETGENPlugin_NETGEN_2D{faceAlgId, &generator};
      mesh->AddHypothesis(compound, faceAlgId);
      hypotheses.emplace_back(alg);
      auto hypId = generator.GetANewId();
      auto hyp = new NETGENPlugin_Hypothesis_2D{hypId, &generator};
      hyp->SetMinSize(minMeshSize);
      hyp->SetMaxSize(maxMeshSize);
      hyp->SetFineness(NETGENPlugin_Hypothesis::UserDefined);
      hyp->SetGrowthRate(0.2);
      hyp->SetNbSegPerEdge(1);
      hyp->SetNbSegPerRadius(1);
      mesh->AddHypothesis(compound, hypId);
#else
      std::cerr << "Cannot generate quad-dominant mesh on general faces "
                   "without Netgen\n";
#endif
    } else {
      hypotheses.emplace_back(new StdMeshers_MEFISTO_2D{faceAlgId, &generator});
      mesh->AddHypothesis(compound, faceAlgId);
    }
  }

  std::unordered_set<TopoDS_Shape, ShapeMapHasher_Hash, ShapeMapHasher_KeyEqual>
      generated_shapes;
  // Loop through shape map (which contains faces and edges)
  for (auto &shape : map_) {
    // Check if shape is not a child (lower dimension) entity
    if (!isChild(shape.first)) {
      static constexpr std::array<TopAbs_ShapeEnum, 3> shapeTypes{
          TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE};
      // Loop through the shape types
      for (auto &shapeType : shapeTypes) {
        for (TopExp_Explorer explorer{shape.first, shapeType}; explorer.More();
             explorer.Next()) {
          auto &subShape =
              explorer.Current();  // The current shape from explorer
          // Search shape map for explored shape
          auto findIter = map_.find(subShape);
          // If it is in the shape map and not in the generated_shapes set
          if (findIter != map_.end() &&
              generated_shapes.find(subShape) == generated_shapes.end()) {
            // If the NucMesh meshing data (params) exist, set algo
            if (auto nmData =
                    dynamic_cast<NucMeshShapeData *>(findIter->second.get())) {
              nmData->setupAlgos(explorer.Current(), generator, *mesh,
                                 hypotheses, minMeshSize, maxMeshSize);
            }
            generated_shapes.emplace(subShape);
          }
        }
      }
    }
  }
  algosTime_.stop();
  meshingTime_.start();
  generator.Compute(*mesh, compound);
  return mesh;
}

}  // namespace NUCMESH
}  // namespace NEM
