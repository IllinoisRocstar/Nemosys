#include "NucMesh/NucMeshShapeData.H"

#include <SMESHDS_GroupBase.hxx>
#include <SMESHDS_GroupOnFilter.hxx>
#include <SMESHDS_GroupOnGeom.hxx>
#include <SMESH_ControlsDef.hxx>
#include <SMESH_Gen.hxx>
#include <SMESH_Group.hxx>
#include <SMESH_HypoFilter.hxx>
#include <SMESH_Mesh.hxx>
#include <StdMeshers_NumberOfSegments.hxx>
#include <StdMeshers_Quadrangle_2D.hxx>
#include <StdMeshers_Regular_1D.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Shape.hxx>
#ifdef HAVE_NGEN
#  include <NETGENPlugin_NETGEN_2D_ONLY.hxx>
#  include <StdMeshers_QuadranglePreference.hxx>
#else
#  include <iostream>
#endif

#include "Mesh/smeshUtils.H"

namespace NEM {
namespace NUCMESH {

void GroupData::setupAlgos(
    const TopoDS_Shape &shape, SMESH_Gen &generator, SMESH_Mesh &mesh,
    std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps) const {
  if (!groupName.empty()) {
    bool addedToGroup = false;
    int oldGroup = -1;
    bool removeOldGroup = false;
    for (auto group : NEM::MSH::containerWrapper(mesh.GetGroups())) {
      auto groupDS = group->GetGroupDS();
      if (groupName == group->GetName() && groupDS &&
          type == group->GetGroupDS()->GetType()) {
        if (auto otherGeomGroup =
                dynamic_cast<SMESHDS_GroupOnGeom *>(groupDS)) {
          auto otherPred = boost::make_shared<SMESH::Controls::BelongToGeom>();
          otherPred->SetGeom(otherGeomGroup->GetShape());
          otherPred->SetType(type);
          auto thisPred = boost::make_shared<SMESH::Controls::BelongToGeom>();
          thisPred->SetGeom(shape);
          thisPred->SetType(type);
          auto disjunction = boost::make_shared<SMESH::Controls::LogicalOR>();
          disjunction->SetPredicate1(otherPred);
          disjunction->SetPredicate2(thisPred);
          mesh.AddGroup(type, groupName.c_str(), -1, TopoDS_Shape{},
                        disjunction);
          addedToGroup = true;
          oldGroup = group->GetID();
          removeOldGroup = true;
          break;
        } else if (auto otherFilterGroup =
                       dynamic_cast<SMESHDS_GroupOnFilter *>(groupDS)) {
          auto otherPred = otherFilterGroup->GetPredicate();
          auto thisPred = boost::make_shared<SMESH::Controls::BelongToGeom>();
          thisPred->SetGeom(shape);
          thisPred->SetType(type);
          auto disjunction = boost::make_shared<SMESH::Controls::LogicalOR>();
          disjunction->SetPredicate1(otherPred);
          disjunction->SetPredicate2(thisPred);
          otherFilterGroup->SetPredicate(disjunction);
          addedToGroup = true;
          break;
        }
      }
    }
    if (!addedToGroup) { mesh.AddGroup(type, groupName.c_str(), -1, shape); }
    if (removeOldGroup && oldGroup >= 0) { mesh.RemoveGroup(oldGroup); }
  }
}

EdgeSegments::EdgeSegments(std::string groupName, int numSegments)
    : CRTPBase(std::move(groupName)), numSegments_(numSegments) {}

void EdgeSegments::setupAlgos(
    const TopoDS_Shape &shape, SMESH_Gen &generator, SMESH_Mesh &mesh,
    std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps) const {
  this->SideSetEdge::setupAlgos(shape, generator, mesh, generatedHyps);
  auto hypId = generator.GetANewId();
  auto hypothesis = new StdMeshers_NumberOfSegments{hypId, &generator};
  hypothesis->SetNumberOfSegments(numSegments_);
  generatedHyps.emplace_back(hypothesis);
  mesh.AddHypothesis(shape, hypId);
  std::list<const SMESHDS_Hypothesis *> algos;
  SMESH_HypoFilter filter(SMESH_HypoFilter::IsAlgo());
  filter.And(SMESH_HypoFilter::IsApplicableTo(shape));
  mesh.GetHypotheses(mesh.GetShapeToMesh(), filter, algos, false);
  int algoId = -1;
  for (auto &hyp : algos) {
    if (auto algo = dynamic_cast<const SMESH_Algo *>(hyp)) {
      // The method is not virtual, and doesn't allow us to modify algo, so
      // const_cast should be safe
      auto allowed_hyps =
          const_cast<SMESH_Algo *>(algo)->GetCompatibleHypothesis();
      if (std::find(allowed_hyps.begin(), allowed_hyps.end(),
                    hypothesis->GetName()) != allowed_hyps.end()) {
        algoId = algo->GetID();
      }
    }
  }
  if (algoId < 0) {
    algoId = generator.GetANewId();
    auto algo = new StdMeshers_Regular_1D{algoId, &generator};
    generatedHyps.emplace_back(algo);
  }
  mesh.AddHypothesis(shape, algoId);
}

void TriMeshSurface::setupAlgos(
    const TopoDS_Shape &shape, SMESH_Gen &generator, SMESH_Mesh &mesh,
    std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps) const {
  this->GroupData::setupAlgos(shape, generator, mesh, generatedHyps);
}

void QuadMeshSurface::setupAlgos(
    const TopoDS_Shape &shape, SMESH_Gen &generator, SMESH_Mesh &mesh,
    std::vector<std::unique_ptr<SMESH_Hypothesis>> &generatedHyps) const {
  this->GroupData::setupAlgos(shape, generator, mesh, generatedHyps);
  int numEdges = 0;
  for (TopExp_Explorer explorer{shape, TopAbs_EDGE}; explorer.More();
       explorer.Next()) {
    ++numEdges;
  }
  if (numEdges == 4) {
    auto algId = generator.GetANewId();
    auto alg = new StdMeshers_Quadrangle_2D{algId, &generator};
    mesh.AddHypothesis(shape, algId);
    generatedHyps.emplace_back(alg);
  } else {
#ifdef HAVE_NGEN
    auto algId = generator.GetANewId();
    auto alg = new NETGENPlugin_NETGEN_2D_ONLY{algId, &generator};
    mesh.AddHypothesis(shape, algId);
    generatedHyps.emplace_back(alg);
    auto hypId = generator.GetANewId();
    auto hyp = new StdMeshers_QuadranglePreference{hypId, &generator};
    mesh.AddHypothesis(shape, hypId);
    generatedHyps.emplace_back(hyp);
#else
    std::cerr << "Cannot generate quad-dominant mesh on general faces "
                 "without Netgen\n";
#endif
  }
}

}  // namespace NUCMESH
}  // namespace NEM