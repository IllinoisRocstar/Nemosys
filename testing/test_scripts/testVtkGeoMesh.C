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
#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif
#include <gtest/gtest.h>

#include <Mesh/vtkGeoMesh.H>

std::string fileNameLegacyUG;
std::string fileNameLegacyUnsupported;
std::string fileNameXMLvtu;
std::string fileNameHex;
std::string fileNameXMLvtuGeo;

TEST(vtkGeoMesh, ConstructorDefault) { NEM::MSH::vtkGeoMesh vgm{}; }

TEST(vtkGeoMesh, ConstructorVtkUnstructuredGrid) {
  NEM::MSH::vtkGeoMesh vgm{vtkSmartPointer<vtkUnstructuredGrid>::New()};
}

TEST(vtkGeoMesh, ConstructorVtkUnstructuredGridWithGeo) {
  NEM::MSH::vtkGeoMesh vgm{vtkSmartPointer<vtkUnstructuredGrid>::New(), "GeoIds"};
}

TEST(vtkGeoMesh, FactoryReadLegacy) {
  NEM::MSH::vtkGeoMesh *vgm = NEM::MSH::vtkGeoMesh::Read(fileNameLegacyUG);
  vgm->Delete();
}

TEST(vtkGeoMeshDeathTest, FactoryReadLegacyUnsupported) {
  EXPECT_DEATH(NEM::MSH::vtkGeoMesh::Read(fileNameLegacyUnsupported),
               "ERROR:.*");
}

TEST(vtkGeoMesh, FactoryReadXML) {
  NEM::MSH::vtkGeoMesh *vgm = NEM::MSH::vtkGeoMesh::Read(fileNameXMLvtu);
  vgm->Delete();
}

//TEST(vtkGeoMeshDeathTest, FactoryReadXMLUnsupported) {
//  EXPECT_DEATH(
//      NEM::MSH::vtkGeoMesh vgm = NEM::MSH::vtkGeoMesh::Read(fileNameXMLvti),
//      "ERROR:.*");
//}

TEST(vtkGeoMeshDeathTest, FactoryReadUnsupportedExtension) {
  EXPECT_DEATH(NEM::MSH::vtkGeoMesh::Read("a.notVTK"), "ERROR:.*notVTK.*");
}

TEST(vtkGeoMesh, FactoryReadXMLWithGeo) {
  NEM::MSH::vtkGeoMesh *vgm = NEM::MSH::vtkGeoMesh::Read(fileNameXMLvtuGeo,
                                                         "material");
  vgm->Delete();
}

TEST(vtkGeoMesh, Write) {
  vtkSmartPointer<vtkUnstructuredGrid> mesh =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  // 0-7: for a hex
  points->InsertNextPoint(0, 0, 0);
  points->InsertNextPoint(1, 0, 0);
  points->InsertNextPoint(1, 1, 0);
  points->InsertNextPoint(0, 1, 0);
  points->InsertNextPoint(0, 0, 1);
  points->InsertNextPoint(1, 0, 1);
  points->InsertNextPoint(1, 1, 1);
  points->InsertNextPoint(0, 1, 1);

  mesh->SetPoints(points);
  vtkIdType hex[8] = {0, 1, 2, 3, 4, 5, 6, 7};
  mesh->InsertNextCell(VTK_HEXAHEDRON, 8, hex);

  NEM::MSH::vtkGeoMesh vgm{mesh};

  vgm.write("hex.vtu");

  NEM::MSH::vtkGeoMesh *vgm2 = NEM::MSH::vtkGeoMesh::Read("hex.vtu");

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  vgm2->getCell(0, cell);

  EXPECT_TRUE(cell->GetCellType() == VTK_HEXAHEDRON);
  for (int i = 0; i != 7; ++i) {
    EXPECT_TRUE(cell->GetPointIds()->GetId(i) == i);
  }

  vgm2->Delete();
}

TEST(vtkGeoMesh, reconstructGeoHex) {
  NEM::MSH::vtkGeoMesh *vgm = NEM::MSH::vtkGeoMesh::Read(fileNameHex);
  vgm->reconstructGeo();
  vgm->Delete();
}

//TEST(geoMeshBase, computeDiscreteGeoFromVTU) {
//  NEM::MSH::vtkGeoMesh *vgm = NEM::MSH::vtkGeoMesh::Read("burned/volumetric_0.vtu");
//
//  vgm->reconstructGeo();
//
//  int fieldTag = gmsh::model::mesh::field::add("Cylinder");
//  gmsh::model::mesh::field::setNumber(fieldTag, "Radius", 0.0175016280264);
//  gmsh::model::mesh::field::setNumber(fieldTag, "VIn", 0.001);
//  gmsh::model::mesh::field::setNumber(fieldTag, "VOut", 0.01);
//  gmsh::model::mesh::field::setNumber(fieldTag, "XAxis", 0.316356986761);
//  gmsh::model::mesh::field::setNumber(fieldTag, "YAxis", 0);
//  gmsh::model::mesh::field::setNumber(fieldTag, "ZAxis", 0);
//  gmsh::model::mesh::field::setNumber(fieldTag, "XCenter", 0);
//  gmsh::model::mesh::field::setNumber(fieldTag, "YCenter", 0);
//  gmsh::model::mesh::field::setNumber(fieldTag, "ZCenter", 0);
//  gmsh::model::mesh::field::setAsBackgroundMesh(fieldTag);
//
//  gmsh::write("burned/volumetric_0_mid.msh");
//  gmsh::model::mesh::generate(3);
//
//  gmsh::write("burned/volumetric_0_remeshed.msh");
//  gmsh::write("burned/volumetric_0_remeshed.vtk");
//
//  vgm = NEM::MSH::vtkGeoMesh::Read("burned/volumetric_0_remeshed.vtk");
//  vgm->write("burned/volumetric_0_remeshed.vtu");
//
//  vgm->Delete();
//}

/*
TEST(geoMeshBase, packFromVTU) {
//  std::string filename = "two_mesh";
//  std::string filename = "packing/TwoSpheres";
//  std::string filename = "cyl2box/cyl2box";
//  std::string filename = "last_packing/Coarser_MixedShapes";
  std::string filename = "nucmesh_sample/nucmesh";

  NEM::MSH::vtkGeoMesh *vgm = NEM::MSH::vtkGeoMesh::Read(filename + ".vtu", "PhysGrpId");

  vgm->reconstructGeo();

  gmsh::write(filename + "_mid.msh");

  int fieldTag = gmsh::model::mesh::field::add("Cylinder");
  gmsh::model::mesh::field::setNumber(fieldTag, "Radius", 8);
  gmsh::model::mesh::field::setNumber(fieldTag, "VIn", 0.1);
  gmsh::model::mesh::field::setNumber(fieldTag, "VOut", 0.2);
  gmsh::model::mesh::field::setNumber(fieldTag, "XAxis", 0);
  gmsh::model::mesh::field::setNumber(fieldTag, "YAxis", 0);
  gmsh::model::mesh::field::setNumber(fieldTag, "ZAxis", 3.1);
  gmsh::model::mesh::field::setNumber(fieldTag, "XCenter", 0);
  gmsh::model::mesh::field::setNumber(fieldTag, "YCenter", 0);
  gmsh::model::mesh::field::setNumber(fieldTag, "ZCenter", 0);
//  gmsh::model::mesh::field::setAsBackgroundMesh(fieldTag);

  int restrictTag = gmsh::model::mesh::field::add("Restrict");
  gmsh::model::mesh::field::setNumber(restrictTag, "IField", fieldTag);

  std::vector<int> threeIds;
  std::vector<double> vertIds, edgeIds, surfIds, regIds;
  gmsh::vectorpair dim3Tags, dim2Tags, dim1Tags, dim0Tags;
  gmsh::model::getEntitiesForPhysicalGroup(3, 4, threeIds);
  for (const auto &threeId : threeIds) dim3Tags.emplace_back(3, threeId);

  gmsh::model::getBoundary(dim3Tags, dim2Tags, false, false, false);
  gmsh::model::getBoundary(dim2Tags, dim1Tags, false, false, false);
  gmsh::model::getBoundary(dim1Tags, dim0Tags, false, false, false);
  for (const auto &dim3Tag : dim3Tags) regIds.emplace_back(dim3Tag.second);
  for (const auto &dim2Tag : dim2Tags) surfIds.emplace_back(dim2Tag.second);
  for (const auto &dim1Tag : dim1Tags) edgeIds.emplace_back(dim1Tag.second);
  for (const auto &dim0Tag : dim0Tags) vertIds.emplace_back(dim0Tag.second);
  gmsh::model::mesh::field::setNumbers(restrictTag, "RegionsList", regIds);
  gmsh::model::mesh::field::setNumbers(restrictTag, "FacesList", surfIds);
  gmsh::model::mesh::field::setNumbers(restrictTag, "EdgesList", edgeIds);
  gmsh::model::mesh::field::setNumbers(restrictTag, "VerticesList", vertIds);
  gmsh::model::mesh::field::setAsBackgroundMesh(restrictTag);

  gmsh::option::setNumber("Mesh.Algorithm", 6);
  gmsh::option::setNumber("Mesh.Algorithm3D", 10);
  gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0);
  gmsh::option::setNumber("Mesh.CharacteristicLengthMax", 0.2);
  gmsh::model::mesh::generate(3);

  gmsh::write(filename + "_remeshed.msh");
  gmsh::write(filename + "_remeshed.vtk");

  vgm = NEM::MSH::vtkGeoMesh::Read(filename + "_remeshed.vtk");
  vgm->write(filename + "_remeshed.vtu");

  vgm->Delete();
}
*/

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  fileNameLegacyUG = argv[1];
  fileNameLegacyUnsupported = argv[2];
  fileNameXMLvtu = argv[3];
  fileNameHex = argv[3];
  fileNameXMLvtuGeo = argv[4];

  int res = RUN_ALL_TESTS();

  return res;
}
