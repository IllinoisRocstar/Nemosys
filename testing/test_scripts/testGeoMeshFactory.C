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
#include <gtest/gtest.h>
#include <Mesh/geoMeshFactory.H>

#include <Mesh/vtkGeoMesh.H>
#include <Mesh/oshGeoMesh.H>
#include <Mesh/exoGeoMesh.H>

#ifdef HAVE_GMSH
#include <Mesh/gmshGeoMesh.H>
#endif

#ifdef HAVE_CFMSH
#  include <Mesh/foamGeoMesh.H>
#endif

// Tests for NEM::MSH::Read accepting files of all geoMeshBase child classes

std::string arg_vtkMesh;
std::string arg_gmshMesh;
std::string arg_oshMesh;
std::string arg_exoMesh;
std::string arg_ofMesh;

TEST(geoMeshFactoryDeathTest, ReadNotSupported) {
  std::string fileName = "fileExt.notSupported";

  EXPECT_DEATH(NEM::MSH::Read(fileName),
               "File extension .notSupported is not supported.");
}

TEST(geoMeshFactory, ReadVTK) {
  NEM::MSH::geoMeshBase *vtk1 = NEM::MSH::Read(arg_vtkMesh);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::vtkGeoMesh *>(vtk1));
  delete vtk1;

  NEM::MSH::geoMeshBase *vtk2 =
      NEM::MSH::Read(arg_vtkMesh, NEM::MSH::MeshType::VTK_GEO_MESH);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::vtkGeoMesh *>(vtk2));
  delete vtk2;
}

#ifdef HAVE_GMSH
TEST(geoMeshFactory, ReadGMSH) {
  NEM::MSH::geoMeshBase *gmsh1 = NEM::MSH::Read(arg_gmshMesh);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::gmshGeoMesh *>(gmsh1));
  delete gmsh1;

  NEM::MSH::geoMeshBase *gmsh2 =
      NEM::MSH::Read(arg_gmshMesh, NEM::MSH::MeshType::GMSH_GEO_MESH);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::gmshGeoMesh *>(gmsh2));
  delete gmsh2;
}
#endif

TEST(geoMeshFactory, ReadOSH) {
  NEM::MSH::geoMeshBase *osh1 = NEM::MSH::Read(arg_oshMesh);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::oshGeoMesh *>(osh1));
  delete osh1;

  NEM::MSH::geoMeshBase *osh2 =
      NEM::MSH::Read(arg_oshMesh, NEM::MSH::MeshType::OSH_GEO_MESH);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::oshGeoMesh *>(osh2));
  delete osh2;
}

TEST(geoMeshFactory, ReadEXO) {
  NEM::MSH::geoMeshBase *exo1 = NEM::MSH::Read(arg_exoMesh);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::exoGeoMesh *>(exo1));
  delete exo1;

  NEM::MSH::geoMeshBase *exo2 =
      NEM::MSH::Read(arg_exoMesh, NEM::MSH::MeshType::EXO_GEO_MESH);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::exoGeoMesh *>(exo2));
  delete exo2;
}
#ifdef HAVE_CFMSH
TEST(geoMeshFactory, ReadFOAM) {
  NEM::MSH::geoMeshBase *ofMesh1 = NEM::MSH::Read(arg_ofMesh);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::foamGeoMesh *>(ofMesh1));
  delete ofMesh1;

  NEM::MSH::geoMeshBase *ofMesh2 =
      NEM::MSH::Read(arg_ofMesh, NEM::MSH::MeshType::FOAM_GEO_MESH);
  EXPECT_TRUE(dynamic_cast<NEM::MSH::foamGeoMesh *>(ofMesh2));
  delete ofMesh2;
}
#endif

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  assert(argc == 6);
  arg_vtkMesh = argv[1];
  arg_gmshMesh = argv[2];
  arg_oshMesh = argv[3];
  arg_exoMesh = argv[4];
  arg_ofMesh = argv[5];

  int res = RUN_ALL_TESTS();

  return res;
}
