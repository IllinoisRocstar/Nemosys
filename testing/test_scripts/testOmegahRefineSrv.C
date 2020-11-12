#include <gtest.h>

#include "omegahRefineSrv.H"

#include <geoMeshFactory.H>
#include <vtkSmartPointer.h>
#include <Omega_h_build.hpp>

#include "diffMesh.H"
#include "oshGeoMesh.H"
#include "vtkGeoMesh.H"

// Test cases for NEM::SRV::omegahRefineSrv

TEST(omegahRefineSrv, ConstructorDefault) { NEM::SRV::omegahRefineSrv sb{}; }

TEST(omegahRefineSrv, ExecuteEmptyMesh) {
  auto sb = vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  auto in = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  auto out = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();

  sb->SetInputDataObject(in);
  sb->Update();
  out = NEM::MSH::oshGeoMesh::SafeDownCast(sb->GetOutputDataObject(0));
}

TEST(omegahRefineSrv, Execute) {
  Omega_h::Mesh boxMesh =
      Omega_h::build_box(NEM::MSH::OmegaHInterface::GetLibrary()->world(),
                         OMEGA_H_SIMPLEX, 1.0, 1.0, 1.0, 10, 10, 10);

  auto sb = vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  auto in = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  in->setOshMesh(&boxMesh);

  sb->SetInputDataObject(in);
  sb->Update();

  auto out = NEM::MSH::oshGeoMesh::SafeDownCast(sb->GetOutput());
  auto gold = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::Read("gold_refined_box.osh"));
  for (auto dim : {0, in->getOshMesh().dim()} ) {
    for (Omega_h::Int i = 0; i < in->getOshMesh().ntags(dim); ++i) {
      auto tag = in->getOshMesh().get_tag(dim, i);
      if (tag->type() == OMEGA_H_REAL) {
        EXPECT_TRUE(out->getOshMesh().has_tag(dim, tag->name()));
      }
    }
  }
  EXPECT_EQ(0, NEM::MSH::diffMesh(gold, out, 1e-9, 0.03));
}

TEST(omegahRefineSrv, RefineUniform) {
  auto in = vtkSmartPointer<NEM::MSH::oshGeoMesh>::Take(
      NEM::MSH::oshGeoMesh::Read("unrefined_beam_tet.osh"));

  auto sb = vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  sb->SetInputDataObject(in);
  sb->Update();

  auto out = NEM::MSH::oshGeoMesh::SafeDownCast(sb->GetOutput());
  auto gold = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::Read("gold_refined_beam_uniform_srv.osh"));
  for (auto dim : {0, in->getOshMesh().dim()} ) {
    for (Omega_h::Int i = 0; i < in->getOshMesh().ntags(dim); ++i) {
      auto tag = in->getOshMesh().get_tag(dim, i);
      if (tag->type() == OMEGA_H_REAL) {
        EXPECT_TRUE(out->getOshMesh().has_tag(dim, tag->name()));
      }
    }
  }
  EXPECT_EQ(0, NEM::MSH::diffMesh(gold, out));
}

// Test that reconstructGeo works for 2D meshes; test adding metric source
TEST(omegahRefineSrv, RefineValueInterface2D) {
  auto mesh = vtkSmartPointer<NEM::MSH::vtkGeoMesh>::Take(
      NEM::MSH::vtkGeoMesh::Read("stage_13_0.vtu", "Mat ID"));
  mesh->reconstructGeo();
  auto in = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  in->takeGeoMesh(mesh);

  auto sb = vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  sb->AddMetricSource(OMEGA_H_VARIATION, 800., "S");
  sb->SetInputDataObject(in);
  sb->Update();

  auto out = NEM::MSH::oshGeoMesh::SafeDownCast(sb->GetOutput());
  auto gold = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::Read("gold_refined_stage_13_0.osh"));
  for (auto dim : {0, in->getOshMesh().dim()} ) {
    for (Omega_h::Int i = 0; i < in->getOshMesh().ntags(dim); ++i) {
      auto tag = in->getOshMesh().get_tag(dim, i);
      if (tag->type() == OMEGA_H_REAL) {
        EXPECT_TRUE(out->getOshMesh().has_tag(dim, tag->name()));
      }
    }
  }
  EXPECT_EQ(0, NEM::MSH::diffMesh(gold, out));
}

// Test that reconstructGeo works for 3D meshes; test specifying transfer
TEST(omegahRefineSrv, RefineValueInterface3D) {
  auto mesh = vtkSmartPointer<NEM::MSH::vtkGeoMesh>::Take(
      NEM::MSH::vtkGeoMesh::Read("stage_9_0.vtu", "Mat ID"));
  mesh->reconstructGeo();
  auto in = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  in->takeGeoMesh(mesh);

  auto sb = vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  sb->AddMetricSource(OMEGA_H_VARIATION, 800., "S");
  sb->AddTransferOpts("S", OMEGA_H_LINEAR_INTERP);
  sb->SetInputDataObject(in);
  sb->Update();

  auto out = NEM::MSH::oshGeoMesh::SafeDownCast(sb->GetOutput());
  auto gold = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::Read("gold_refined_stage_9_0.osh"));
  EXPECT_EQ(0, NEM::MSH::diffMesh(gold, out, 1e-9, 1e-5, 0.02, 0.02));
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  int res = RUN_ALL_TESTS();

  return res;
}
