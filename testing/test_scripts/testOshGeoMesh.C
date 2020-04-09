#include <gtest.h>

#include "Omega_h_build.hpp"
#include "Omega_h_file.hpp"
#include "vtkXMLDataSetWriter.h"
#include "vtkXMLGenericDataObjectReader.h"

#include "oshGeoMesh.H"

TEST(oshGeoMesh, ConstructorDefault) { NEM::MSH::oshGeoMesh ogm{}; }

TEST(oshGeoMesh, ConstructorOmegaHMesh) {
  Omega_h::Mesh oshMesh{};
  NEM::MSH::oshGeoMesh ogm{&oshMesh};
}

TEST(oshGeoMesh, ConstructorOmegaHMeshValid) {
  Omega_h::Library oshLibrary{};
  Omega_h::Mesh oshMesh{&oshLibrary};
  Omega_h::build_box_internal(&oshMesh, OMEGA_H_HYPERCUBE, 1.0, 1.0, 1.0, 2, 2,
                              2);
  NEM::MSH::oshGeoMesh ogm{&oshMesh};
}

TEST(oshGeoMesh, ConstructorOmegaHMeshAndLib) {
  Omega_h::Library oshLibrary{};
  Omega_h::Mesh oshMesh{&oshLibrary};
  NEM::MSH::oshGeoMesh ogm{&oshMesh, &oshLibrary};
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
