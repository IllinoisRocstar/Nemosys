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
#include "Services/NucMeshSrv.H"

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <SMESH_ControlsDef.hxx>
#include <SMESH_Group.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESH_MeshEditor.hxx>
#include <SMESH_TypeDefs.hxx>
#include <SMESHDS_GroupBase.hxx>
#include <SMESHDS_Mesh.hxx>

#include "Mesh/smeshGeoMesh.H"
#include "Mesh/smeshUtils.H"
#include "NucMesh/NucMeshGeo.H"
#include "NucMesh/ShapeBase.H"

namespace NEM {
namespace SRV {

vtkStandardNewMacro(NucMeshSrv)

NucMeshSrv::NucMeshSrv() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

void NucMeshSrv::SetConfiguration(const NucMeshConf &configuration) {
  conf_ = configuration;
  this->Modified();
}

const NucMeshConf &NucMeshSrv::GetConfiguration() const { return conf_; }

int NucMeshSrv::FillOutputPortInformation(int port, vtkInformation *info) {
  if (port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "smeshGeoMesh");
    return 1;
  } else {
    return 0;
  }
}

int NucMeshSrv::RequestData(vtkInformation *request,
                            vtkInformationVector **inputVector,
                            vtkInformationVector *outputVector) {
  NEM::NUCMESH::NucMeshGeo geo(0);
  for (auto &pattern : conf_.geometryAndMesh) {
    if (pattern) {
      NEM::NUCMESH::ShapeBase::mergeGeo(geo, pattern->createGeo());
    }
  }
  auto mesh = geo.computeMesh(*this->conf_.generator);
  std::array<TIDSortedElemSet, 2> elems_nodes;
  {
    auto elemContainer =
        NEM::MSH::containerWrapper(mesh->GetMeshDS()->elementsIterator());
    std::copy(elemContainer.begin(), elemContainer.end(),
              std::inserter(elems_nodes[0], elems_nodes[0].end()));
    auto nodesContainer =
        NEM::MSH::containerWrapper(mesh->GetMeshDS()->nodesIterator());
    std::copy(nodesContainer.begin(), nodesContainer.end(),
              std::inserter(elems_nodes[1], elems_nodes[1].end()));
  }
  if (!conf_.extrudeSteps.empty()) {
    std::set<int> groupsToRemove;
    // Remove all old face groups
    for (auto group : NEM::MSH::containerWrapper(mesh->GetGroups())) {
      if (group->GetGroupDS()->GetType() == SMDSAbs_Face) {
        groupsToRemove.emplace(group->GetID());
      }
    }
    const std::string bottomName = "Bottom";
    mesh->AddGroup(SMDSAbs_Face, bottomName.c_str(), -1,
                   mesh->GetShapeToMesh());
    Handle(TColStd_HSequenceOfReal)
        extrudeStepsSeq{new TColStd_HSequenceOfReal{}};
    for (const auto &step : conf_.extrudeSteps) {
      extrudeStepsSeq->Append(step);
    }
    SMESH_MeshEditor::ExtrusParam extrusParam{
        gp_Dir{gp_XYZ(0, 0, 1)}, extrudeStepsSeq,
        SMESH_MeshEditor::EXTRUSION_FLAG_BOUNDARY |
            SMESH_MeshEditor::EXTRUSION_FLAG_GROUPS};
    SMESH_MeshEditor::TTElemOfElemListMap newElems;
    SMESH_MeshEditor editor{mesh.get()};
    auto newGroups =
        editor.ExtrusionSweep(elems_nodes.data(), extrusParam, newElems);
    for (auto &newGroupID : *newGroups) {
      auto group = mesh->GetGroup(newGroupID);
      auto type = group->GetGroupDS()->GetType();
      if (type == SMDSAbs_Volume) {
        // Remove "Bottom_extruded"
        const std::string tempAllExtruded = bottomName + "_extruded";
        if (tempAllExtruded.compare(0, tempAllExtruded.size(),
                                    group->GetName()) == 0) {
          groupsToRemove.emplace(group->GetID());
        }
      } else if (type == SMDSAbs_Face) {
        const std::string tempAllTop = bottomName + "_top";
        if (tempAllTop.compare(0, tempAllTop.size(), group->GetName()) == 0) {
          group->SetName("Top");
        } else {
          std::string groupName = group->GetName();
          const std::string suffix = "_top";
          // If ends in _top, remove it (but keep the extruded side sets)
          if (groupName.size() > suffix.size() &&
              groupName.compare(groupName.size() - suffix.size(), suffix.size(),
                                suffix) == 0) {
            groupsToRemove.emplace(group->GetID());
          }
        }
      }
    }
    for (auto &groupID : groupsToRemove) { mesh->RemoveGroup(groupID); }
  }

  MSH::smeshGeoMesh *output = MSH::smeshGeoMesh::SafeDownCast(
      outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  output->setSMeshMesh(std::move(mesh), std::move(conf_.generator));
  return 1;
}

}  // namespace SRV
}  // namespace NEM