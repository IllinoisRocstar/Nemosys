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
#include "Mesh/smeshGeoMesh.H"

#include <algorithm>
#include <array>
#include <cassert>

#include <SMESH_Gen.hxx>
#include <SMESH_Group.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESHDS_Group.hxx>
#include <SMESHDS_Mesh.hxx>

#include <vtkCellIterator.h>
#include <vtkCellType.h>
#include <vtkDataSet.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>

#include "Mesh/smeshUtils.H"

namespace {

void insertWithGroup(const std::vector<SMDS_MeshNode *> &smdsNodes,
                     SMESHDS_Mesh *meshDS, SMESH_Mesh &outMesh,
                     vtkDataSet *dataset, vtkDataArray *entArr,
                     vtkStringArray *groupNames) {
  std::vector<vtkIdType> pointIdBuffer;
  std::map<int, SMESH_Group *> groups;
  for (auto cellIter =
           vtkSmartPointer<vtkCellIterator>::Take(dataset->NewCellIterator());
       !cellIter->IsDoneWithTraversal(); cellIter->GoToNextCell()) {
    pointIdBuffer.clear();
    auto vtkPointIds = cellIter->GetPointIds();
    auto cellDim = cellIter->GetCellDimension();
    SMDS_MeshCell *outCell;
    if (cellDim == 0) {
      outCell = meshDS->Add0DElement(smdsNodes[vtkPointIds->GetId(0)]);
    } else if (cellDim == 1) {
      auto cellType = cellIter->GetCellType();
      if (cellType == VTK_LINE) {
        outCell = meshDS->AddEdge(smdsNodes[vtkPointIds->GetId(0)],
                                  smdsNodes[vtkPointIds->GetId(1)]);
      } else if (cellType == VTK_QUADRATIC_EDGE) {
        outCell = meshDS->AddEdge(smdsNodes[vtkPointIds->GetId(0)],
                                  smdsNodes[vtkPointIds->GetId(1)],
                                  smdsNodes[vtkPointIds->GetId(2)]);
      }
    } else if (cellDim == 2 || cellDim == 3) {
      for (int i = 0; i < vtkPointIds->GetNumberOfIds(); ++i) {
        pointIdBuffer.emplace_back(
            smdsNodes[vtkPointIds->GetId(i)]->GetVtkID());
      }
      if (cellDim == 2) {
        outCell = meshDS->AddFaceFromVtkIds(pointIdBuffer);
      } else {
        assert(cellDim == 3);
        outCell = meshDS->AddVolumeFromVtkIds(pointIdBuffer);
      }
    }
    if (!outCell) {
      std::cerr << "Failed to add cell with vtk type "
                << cellIter->GetCellType() << '\n';
    }
    if (entArr) {
      auto gmGroupId =
          static_cast<int>(entArr->GetComponent(cellIter->GetCellId(), 0));
      auto &group = groups[gmGroupId];
      if (!group) {
        std::string groupName;
        if (groupNames && gmGroupId >= 0 &&
            gmGroupId < groupNames->GetNumberOfValues()) {
          groupName = groupNames->GetValue(gmGroupId);
        }
        if (groupName.empty()) {
          groupName = "Group " + std::to_string(gmGroupId);
        }
        group = outMesh.AddGroup(SMDSAbs_All, groupName.c_str());
      }
      if (auto manualGroup =
              dynamic_cast<SMESHDS_Group *>(group->GetGroupDS())) {
        manualGroup->Add(outCell);
      }
    }
  }
}

}  // namespace

namespace NEM {
namespace MSH {

vtkStandardNewMacro(smeshGeoMesh)

smeshGeoMesh::smeshGeoMesh()
    : geoMeshBase(), gen_{new SMESH_Gen{}}, mesh_(gen_->CreateMesh(false)) {}

// Not in header so we can forward declare SMESH_Mesh and SMESH_Gen
smeshGeoMesh::~smeshGeoMesh() = default;

void smeshGeoMesh::write(const std::string &fileName) {
  std::cerr << "Currently unsupported.\n";
}

void smeshGeoMesh::report(std::ostream &out) const { geoMeshBase::report(out); }

void smeshGeoMesh::resetNative() {
  // Careful! Delete mesh before deleting gen!
  mesh_.reset();
  gen_.reset(new SMESH_Gen{});
  mesh_.reset(gen_->CreateMesh(false));
  smeshGeoMesh::GMToSMESH(getGeoMesh(), *mesh_);
}

void smeshGeoMesh::setSMeshMesh(std::unique_ptr<SMESH_Mesh> &&mesh,
                                std::shared_ptr<SMESH_Gen> gen) {
  mesh_ = std::move(mesh);
  gen_ = std::move(gen);
  setGeoMesh(SmeshToGM(*mesh_));
}

const SMESH_Mesh &smeshGeoMesh::getSMESHMesh() const { return *mesh_; }

void smeshGeoMesh::GMToSMESH(const GeoMesh &geoMesh, SMESH_Mesh &outMesh) {
  auto grid = geoMesh.mesh;
  if (grid->GetNumberOfCells() == 0) { return; }
  auto linkArr = grid->GetCellData()->GetArray(geoMesh.link.c_str());
  auto points = grid->GetPoints();
  auto meshDS = outMesh.GetMeshDS();
  std::vector<SMDS_MeshNode *> smdsNodes;  // indexed by geoMesh id
  smdsNodes.reserve(points->GetNumberOfPoints());
  {
    std::array<double, 3> pointBuffer{};
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
      grid->GetPoint(i, pointBuffer.data());
      smdsNodes.emplace_back(
          meshDS->AddNode(pointBuffer[0], pointBuffer[1], pointBuffer[2]));
    }
  }
  insertWithGroup(smdsNodes, meshDS, outMesh, grid, linkArr,
                  vtkStringArray::SafeDownCast(
                      grid->GetFieldData()->GetArray("Group Names")));
  insertWithGroup(smdsNodes, meshDS, outMesh, geoMesh.sideSet.sides,
                  geoMesh.sideSet.getGeoEntArr(),
                  geoMesh.sideSet.getSideSetNames());
}

// Note not hardened to support cells belonging to multiple groups
geoMeshBase::GeoMesh smeshGeoMesh::SmeshToGM(SMESH_Mesh &mesh) {
  vtkNew<vtkUnstructuredGrid> gmGrid;
  vtkNew<vtkPolyData> sideSetSides;
  int numGMCells;
  int maxDim;
  if ((numGMCells = mesh.NbVolumes()) > 0) {
    maxDim = 3;
  } else if ((numGMCells = mesh.NbFaces()) > 0) {
    maxDim = 2;
  } else if ((numGMCells = mesh.NbEdges()) > 0) {
    maxDim = 1;
  } else {
    return {gmGrid, {}, {}, {}};
  }
  auto smeshGrid = mesh.GetMeshDS()->GetGrid();
  gmGrid->Allocate(numGMCells);
  gmGrid->SetPoints(smeshGrid->GetPoints());
  sideSetSides->Allocate();
  sideSetSides->SetPoints(smeshGrid->GetPoints());
  vtkNew<vtkStringArray> groupNamesGrid;
  groupNamesGrid->SetName("Group Names");
  vtkNew<vtkStringArray> groupNamesSides;
  // indexed by vtkID on mesh.GetMeshDS()->GetGrid()
  std::vector<int> groups(smeshGrid->GetNumberOfCells(), -1);
  for (auto group : NEM::MSH::containerWrapper(mesh.GetGroups())) {
    auto groupId = group->GetID();
    if (groupId >= 0) {
      for (auto elem :
           NEM::MSH::containerWrapper(group->GetGroupDS()->GetElements())) {
        auto vtkId = elem->GetVtkID();
        if (vtkId >= 0) {
          auto &cellGroup = groups[vtkId];
          if (cellGroup >= 0) {
            // For now, the behavior is to choose the "latest" group.
            std::cerr << "Warning! SMESH Element found on multiple groups. "
                         "Behavior is undefined.\n";
          }
          cellGroup = groupId;
        }
      }
      auto elemType = group->GetGroupDS()->GetType();
      auto groupDim =
          elemType == SMDSAbs_Volume                                    ? 3
          : elemType == SMDSAbs_Face                                    ? 2
          : (elemType == SMDSAbs_Node || elemType == SMDSAbs_0DElement) ? 1
                                                                        : -1;
      if (groupDim == maxDim || groupDim == -1) {
        groupNamesGrid->InsertValue(groupId, group->GetName());
      }
      if (groupDim == maxDim - 1 || groupDim == -1) {
        groupNamesSides->InsertValue(groupId, group->GetName());
      }
    }
  }
  vtkNew<vtkIntArray> groupsArr;
  groupsArr->SetName("GeoEnt");
  groupsArr->SetNumberOfComponents(1);
  groupsArr->SetNumberOfTuples(numGMCells);
  groupsArr->FillComponent(0, -1);
  vtkNew<vtkIntArray> sideSetEntArr;
  bool geoMeshHasLink;
  for (auto cellIter =
           vtkSmartPointer<vtkCellIterator>::Take(smeshGrid->NewCellIterator());
       !cellIter->IsDoneWithTraversal(); cellIter->GoToNextCell()) {
    auto cellDim = cellIter->GetCellDimension();
    if (cellDim == maxDim) {
      auto gmId = gmGrid->InsertNextCell(cellIter->GetCellType(),
                                         cellIter->GetPointIds());
      auto group = groups[cellIter->GetCellId()];
      if (group >= 0) {
        geoMeshHasLink = true;
        groupsArr->SetComponent(gmId, 0, group);
      }
    } else if (cellDim == maxDim - 1) {
      auto group = groups[cellIter->GetCellId()];
      if (group >= 0) {
        sideSetSides->InsertNextCell(cellIter->GetCellType(),
                                     cellIter->GetPointIds());
        sideSetEntArr->InsertNextValue(group);
      }
    }
  }
  std::string link;
  if (geoMeshHasLink) {
    gmGrid->GetCellData()->AddArray(groupsArr);
    gmGrid->GetFieldData()->AddArray(groupNamesGrid);
    link = groupsArr->GetName();
  }
  return {gmGrid,
          {},
          std::move(link),
          sideSetSides->GetNumberOfCells() > 0
              ? SideSet(sideSetSides, sideSetEntArr, nullptr, nullptr,
                        groupNamesSides)
              : SideSet{}};
}

}  // namespace MSH
}  // namespace NEM