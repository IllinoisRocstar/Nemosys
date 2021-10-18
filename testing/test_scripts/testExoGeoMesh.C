#include <gtest/gtest.h>
#include <vtkInformation.h>
#include <vtkInformationStringKey.h>

#include <Mesh/diffMesh.H>
#include <Mesh/exoGeoMesh.H>
#include <Mesh/geoMeshFactory.H>

// Test cases for the NEM::MSH::exoGeoMesh class, which handles Exodus II files

std::string arg_fName1;
std::string arg_fName2;
std::string arg_fName3;
std::string arg_fName4;
std::string arg_fName5;

TEST(exoGeoMeshDeathTest, ReadFakeFile) {
  EXPECT_DEATH(vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
                   NEM::MSH::exoGeoMesh::Read("fake_file")),
               "Malformed output of vtkExodusIIReader.");
}

TEST(exoGeoMesh, ConstructorDefault) { NEM::MSH::exoGeoMesh egm{}; }

// Test exoGeoMesh::stitch, also test large number of side sets
TEST(exoGeoMesh, Stitch) {
  auto mesh1 = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(arg_fName1));
  auto mesh1Cells = mesh1->getNumberOfCells();
  auto mesh1Points = mesh1->getNumberOfPoints();
  std::string physGroupName = "PhysicalGroup";
  mesh1->addElemBlockProperty(physGroupName);
  for (const auto &blockId : mesh1->getElemBlockIds()) {
    mesh1->setElemBlockProperty(physGroupName, blockId, 1);
  }
  mesh1->setPhysGrpPropertyName(physGroupName);
  EXPECT_EQ(physGroupName, mesh1->getPhysGrpPropertyName());
  mesh1->reconstructGeo();
  int totalSideSetSides = 0;
  for (const auto &id : mesh1->getSideSetIds()) {
    totalSideSetSides += mesh1->getSideSet(id).first.size();
  }
  EXPECT_EQ(13320, totalSideSetSides);
  auto mesh2 = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(arg_fName2));
  auto mesh2Cells = mesh2->getNumberOfCells();
  auto mesh2Points = mesh2->getNumberOfPoints();
  mesh2->reconstructGeo();
  mesh1->stitch(*mesh2);
  mesh1->reconstructGeo();
  EXPECT_EQ(mesh1->getNumberOfCells(), mesh1Cells + mesh2Cells);
  EXPECT_LE(mesh1->getNumberOfPoints(), mesh1Points + mesh2Points);
  EXPECT_GT(mesh1->getNumberOfPoints(), mesh1Points);
  EXPECT_GT(mesh1->getNumberOfPoints(), mesh2Points);
  mesh1->write("stitched.exo");
  auto mesh3 = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read("stitched.exo"));
  EXPECT_EQ(0, NEM::MSH::diffMesh(mesh1, mesh3));
}

TEST(exoGeoMesh, ReadFile) {
  auto mesh = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(arg_fName3));
  EXPECT_EQ("SimpleTitle", mesh->getTitle());
  auto elemBlockIds = mesh->getElemBlockIds();
  EXPECT_EQ(2, elemBlockIds.size());
  EXPECT_EQ(elemBlockIds.size(), mesh->getNumElemBlocks());
  auto properties = mesh->getElemBlockPropertyNames();
  EXPECT_EQ(properties.size(), 2);
  std::set<std::string> truePropertyNames{"SameProp", "DifferentProp"};
  for (const auto &id : elemBlockIds) {
    auto name = mesh->getElemBlockName(id);
    EXPECT_EQ(name, mesh->getElemBlockNames().at(id));
    auto cells = mesh->getElemBlock(id);
    if (name == "FirstBlock") {
      EXPECT_EQ(cells.size(), 2);
      for (const auto &cellId : cells) {
        EXPECT_EQ(VTK_HEXAHEDRON, mesh->getCell(cellId)->GetCellType());
      }
      EXPECT_EQ(16, mesh->getElemBlockProperty("DifferentProp", id));
    } else {
      EXPECT_EQ("SecondBlock", name);
      EXPECT_EQ(cells.size(), 1);
      for (const auto &cellId : cells) {
        EXPECT_EQ(VTK_WEDGE, mesh->getCell(cellId)->GetCellType());
      }
      EXPECT_EQ(17, mesh->getElemBlockProperty("DifferentProp", id));
    }
    for (const auto &cellId : cells) {
      EXPECT_EQ(id, mesh->getElemBlockId(cellId));
    }
    EXPECT_EQ(15, mesh->getElemBlockProperty("SameProp", id));
  }
  auto sideSetIds = mesh->getSideSetIds();
  EXPECT_EQ(1, sideSetIds.size());
  EXPECT_EQ(sideSetIds.size(), mesh->getNumSideSets());
  EXPECT_EQ("UniqueSideSet", mesh->getSideSetName(sideSetIds[0]));
  EXPECT_EQ(mesh->getSideSetName(sideSetIds[0]),
            mesh->getSideSetNames().at(sideSetIds[0]));
  auto sideSet = mesh->getSideSet(sideSetIds[0]);
  std::set<std::pair<int, int>> setSides;
  for (std::size_t i = 0; i < sideSet.first.size(); ++i) {
    setSides.emplace(sideSet.first[i], sideSet.second[i]);
  }
  auto trueSetSides =
      std::set<std::pair<int, int>>{{0, 4}, {0, 0}, {1, 0}, {2, 4}};
  EXPECT_EQ(trueSetSides, setSides);
  auto nodeSetIds = mesh->getNodeSetIds();
  EXPECT_EQ(1, nodeSetIds.size());
  EXPECT_EQ(nodeSetIds.size(), mesh->getNumNodeSets());
  auto nodeSet = mesh->getNodeSet(nodeSetIds[0]);
  EXPECT_EQ("OnlyNodeSet", mesh->getNodeSetName(nodeSetIds[0]));
  EXPECT_EQ(mesh->getNodeSetName(nodeSetIds[0]),
            mesh->getNodeSetNames().at(nodeSetIds[0]));
  auto trueSetNodes = std::set<int>{0, 1, 2, 3, 9, 10, 12, 13};
  EXPECT_EQ(trueSetNodes, std::set<int>(nodeSet.begin(), nodeSet.end()));
}

// Test exoGeoMesh::write
TEST(exoGeoMesh, ReadCopyFile) {
  auto mesh = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(arg_fName3));
  mesh->write("copy" + arg_fName3);
  auto meshCopy = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read("copy" + arg_fName3));
  EXPECT_EQ(0, NEM::MSH::diffMesh(mesh, meshCopy));
  EXPECT_EQ(mesh->getTitle(), meshCopy->getTitle());
  auto sideSetIds = mesh->getSideSetIds();
  EXPECT_EQ(sideSetIds, meshCopy->getSideSetIds());
  for (const auto &id : sideSetIds) {
    EXPECT_EQ(mesh->getSideSetName(id), meshCopy->getSideSetName(id));
    EXPECT_EQ(mesh->getSideSet(id), meshCopy->getSideSet(id));
  }
  auto elemBlockIds = mesh->getElemBlockIds();
  EXPECT_EQ(elemBlockIds, meshCopy->getElemBlockIds());
  auto properties = mesh->getElemBlockPropertyNames();
  EXPECT_EQ(properties.size(), meshCopy->getElemBlockPropertyNames().size());
  for (const auto &id : elemBlockIds) {
    EXPECT_EQ(mesh->getElemBlockName(id), meshCopy->getElemBlockName(id));
    EXPECT_EQ(mesh->getElemBlock(id), meshCopy->getElemBlock(id));
    for (const auto &property : properties) {
      EXPECT_EQ(mesh->getElemBlockProperty(property, id),
                meshCopy->getElemBlockProperty(property, id));
    }
  }
  auto nodeSetIds = mesh->getNodeSetIds();
  EXPECT_EQ(nodeSetIds, meshCopy->getNodeSetIds());
  for (const auto &id : nodeSetIds) {
    EXPECT_EQ(mesh->getNodeSetName(id), meshCopy->getNodeSetName(id));
    EXPECT_EQ(mesh->getNodeSet(id), meshCopy->getNodeSet(id));
  }
}

// Test getters and setters
TEST(exoGeoMesh, CopyByParts) {
  auto mesh = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(arg_fName3));
  auto meshCopy = vtkSmartPointer<NEM::MSH::exoGeoMesh>::New();
  meshCopy->setTitle(mesh->getTitle());
  EXPECT_EQ(mesh->getTitle(), meshCopy->getTitle());
  auto vtkMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkMesh->Allocate(mesh->getNumberOfCells());
  auto vtkMeshPoints = vtkSmartPointer<vtkPoints>::New();
  vtkMeshPoints->Allocate(mesh->getNumberOfPoints());
  std::array<double, 3> coords{};
  for (vtkIdType i = 0; i < mesh->getNumberOfPoints(); ++i) {
    mesh->getPoint(i, &coords);
    vtkMeshPoints->InsertNextPoint(coords.data());
  }
  auto elemBlockIds = mesh->getElemBlockIds();
  for (const auto &id : elemBlockIds) {
    vtkMesh->Reset();
    vtkMesh->SetPoints(vtkMeshPoints);
    auto elemBlock = mesh->getElemBlock(id);
    for (const auto &elem : elemBlock) {
      auto origCell = mesh->getCell(elem);
      auto idList = vtkSmartPointer<vtkIdList>::New();
      idList->Allocate(origCell->GetNumberOfPoints());
      for (int i = 0; i < origCell->GetNumberOfPoints(); ++i) {
        idList->InsertNextId(origCell->GetPointId(i));
      }
      vtkMesh->InsertNextCell(origCell->GetCellType(), idList);
    }
    auto idCopy = meshCopy->addElemBlock(vtkMesh, mesh->getElemBlockName(id));
    auto elemBlockCopy = meshCopy->getElemBlock(idCopy);
    EXPECT_EQ(elemBlock.size(), elemBlockCopy.size());
    for (std::size_t i = 0; i < elemBlock.size(); ++i) {
      auto cell = mesh->getCell(elemBlock[i]);
      auto cellCopy = meshCopy->getCell(elemBlockCopy[i]);
      EXPECT_EQ(cell->GetCellType(), cellCopy->GetCellType());
      for (int j = 0; j < cell->GetNumberOfPoints(); ++j) {
        EXPECT_EQ(cell->GetPointId(j), cellCopy->GetPointId(j));
      }
    }
  }
  EXPECT_EQ(elemBlockIds.size(), meshCopy->getElemBlockIds().size());
  auto sideSetIds = mesh->getSideSetIds();
  for (const auto &id : sideSetIds) {
    auto sideSet = mesh->getSideSet(id);
    auto idCopy = meshCopy->addSideSet(sideSet.first, sideSet.second,
                                       mesh->getSideSetName(id));
    EXPECT_EQ(sideSet, meshCopy->getSideSet(idCopy));
  }
  EXPECT_EQ(sideSetIds.size(), meshCopy->getSideSetIds().size());
  auto nodeSetIds = mesh->getNodeSetIds();
  for (const auto &id : nodeSetIds) {
    auto nodeSet = mesh->getNodeSet(id);
    auto idCopy = meshCopy->addNodeSet(nodeSet, mesh->getNodeSetName(id));
    auto nodeSetCopy = meshCopy->getNodeSet(idCopy);
    EXPECT_EQ(nodeSet.size(), nodeSetCopy.size());
    for (std::size_t i = 0; i < nodeSet.size(); ++i) {
      EXPECT_EQ(mesh->getPoint(nodeSet[i]), meshCopy->getPoint(nodeSetCopy[i]));
    }
  }
  EXPECT_EQ(nodeSetIds.size(), meshCopy->getNodeSetIds().size());
}

// Test I/O of EX_NODAL variables (aka point data) with 1, 3, or 6 components
TEST(exoGeoMesh, ReadWriteNodeData) {
  auto mesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(arg_fName4));
  auto exoGM = vtkSmartPointer<NEM::MSH::exoGeoMesh>::New();
  exoGM->takeGeoMesh(mesh);
  auto exoFileName =
      arg_fName4.substr(0, arg_fName4.find_last_of('.')) + ".exo";
  exoGM->write(exoFileName);
  auto meshCopy = vtkSmartPointer<NEM::MSH::exoGeoMesh>::Take(
      NEM::MSH::exoGeoMesh::Read(exoFileName));
  // Checks EX_NODAL variables (vtkPointData in the vtk representation)
  EXPECT_EQ(0, NEM::MSH::diffMesh(exoGM, meshCopy));
}

#ifdef HAVE_GMSH
TEST(exoGeoMesh, TakeGeoMesh) {
  auto mesh =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(arg_fName5));
  mesh->reconstructGeo();
  auto exoGM = vtkSmartPointer<NEM::MSH::exoGeoMesh>::New();
  exoGM->takeGeoMesh(mesh);
  EXPECT_EQ(6, exoGM->getSideSetIds().size());
  EXPECT_FALSE(exoGM->getPhysGrpPropertyName().empty());
}
#endif

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 6);
  arg_fName1 = argv[1];
  arg_fName2 = argv[2];
  arg_fName3 = argv[3];
  arg_fName4 = argv[4];
  arg_fName5 = argv[5];

  int res = RUN_ALL_TESTS();

  return res;
}
