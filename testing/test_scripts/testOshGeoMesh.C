#include <gtest/gtest.h>
#include "Omega_h_build.hpp"

#include <Mesh/diffMesh.H>
#include <Mesh/oshGeoMesh.H>
#include <Mesh/vtkGeoMesh.H>

// Test cases for NEM::MSH::oshGeoMesh

TEST(oshGeoMesh, ConstructorDefault) { NEM::MSH::oshGeoMesh ogm{}; }

TEST(oshGeoMesh, ConstructorOmegaHMesh) {
  Omega_h::Mesh oshMesh{};
  NEM::MSH::oshGeoMesh ogm{&oshMesh};
}

TEST(oshGeoMesh, ConstructorOmegaHMeshValid) {
  auto oshMesh =
      Omega_h::build_box(NEM::MSH::OmegaHInterface::GetLibrary()->world(),
                         OMEGA_H_HYPERCUBE, 1.0, 1.0, 1.0, 2, 2, 2);
  NEM::MSH::oshGeoMesh ogm{&oshMesh};
}

TEST(oshGeoMesh, ConstructorOmegaHMeshAndLib) {
  Omega_h::Mesh oshMesh{NEM::MSH::OmegaHInterface::GetLibrary().get()};
  NEM::MSH::oshGeoMesh ogm{&oshMesh};
}

TEST(oshGeoMesh, ReadOsh) {
  auto mesh = vtkSmartPointer<NEM::MSH::oshGeoMesh>::Take(
      NEM::MSH::oshGeoMesh::Read("box0.osh"));
  EXPECT_EQ(8, mesh->getOshMesh().nelems());
}

TEST(oshGeoMesh, WriteOsh) {
  auto oshMesh =
      Omega_h::build_box(NEM::MSH::OmegaHInterface::GetLibrary()->world(),
                         OMEGA_H_HYPERCUBE, 1.0, 1.0, 1.0, 2, 2, 2);
  NEM::MSH::oshGeoMesh ogm{&oshMesh};
  ogm.write("boxTestWrite.osh");
  EXPECT_EQ(0, NEM::MSH::diffMesh(
                   &ogm, vtkSmartPointer<NEM::MSH::oshGeoMesh>::Take(
                             NEM::MSH::oshGeoMesh::Read("boxTestWrite.osh"))));
}

// Test GM2osh by calling takeGeoMesh test osh2GM by calling setOshMesh
TEST(oshGeoMesh, GM2osh2GM) {
  auto mesh = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::vtkGeoMesh::Read("two_mesh.vtu"));
  mesh->reconstructGeo();
  auto oshGM1 = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  oshGM1->takeGeoMesh(mesh);
  auto oshGM2 = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  auto oshMeshCopy = oshGM1->getOshMesh();
  oshGM2->setOshMesh(&oshMeshCopy);
  EXPECT_EQ(0, NEM::MSH::diffMesh(oshGM1, oshGM2));
}

/*
TEST(oshGeoMesh, osh2vtk) {
  Omega_h::Library lib{};
  vtkSmartPointer<vtkUnstructuredGrid> ug =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  Omega_h::Mesh osh = Omega_h::gmsh::read("cube41.msh", lib.world());
  ug.TakeReference(NEM::MSH::oshGeoMesh::osh2vtk(&osh));

  vtkSmartPointer<vtkXMLDataSetWriter> writer =
      vtkSmartPointer<vtkXMLDataSetWriter>::New();
  writer->SetFileName("cube41.vtu");
  writer->SetInputData(ug);

  writer->Write();
}
TEST(oshGeoMesh, vtk2osh) {
  Omega_h::Library lib{};
  vtkSmartPointer<vtkUnstructuredGrid> ug =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkXMLGenericDataObjectReader> reader =
      vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
  reader->SetFileName("cube41.vtu");
  reader->Update();

  ug = vtkUnstructuredGrid::SafeDownCast(reader->GetOutput());

  Omega_h::Mesh *osh = NEM::MSH::oshGeoMesh::vtk2osh(ug, &lib);

  Omega_h::binary::write("cube41.osh", osh);
  Omega_h::vtk::write_vtu("cube41_out.vtu", osh);

  delete osh;
}
*/

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  int res = RUN_ALL_TESTS();

  return res;
}
