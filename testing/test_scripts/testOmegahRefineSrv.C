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
#include <Services/omegahRefineSrv.H>

#include <Mesh/geoMeshFactory.H>
#include <vtkSmartPointer.h>
#include <vtkExecutive.h>
#include <vtkCommand.h>
#include <Omega_h_build.hpp>

#include <Mesh/diffMesh.H>
#include <Mesh/oshGeoMesh.H>
#include <Mesh/vtkGeoMesh.H>

// Test cases for NEM::SRV::omegahRefineSrv

class ErrorObserver : public vtkCommand {
 public:
  static ErrorObserver *New() { return new ErrorObserver; }
  void Execute(vtkObject *caller, unsigned long event,
               void *calldata) override {
    if (event == vtkCommand::ErrorEvent) {
      error = true;
      errorMsg = static_cast<const char *>(calldata);
    }
  }
  bool error{false};
  std::string errorMsg{};
};

TEST(omegahRefineSrv, ExecuteEmptyMesh) {
  auto sb = vtkSmartPointer<NEM::SRV::omegahRefineSrv>::New();
  auto in = vtkSmartPointer<NEM::MSH::oshGeoMesh>::New();
  vtkNew<ErrorObserver> observer;

  sb->SetInputDataObject(in);
  sb->GetExecutive()->AddObserver(vtkCommand::AnyEvent, observer);
  sb->Update();
  EXPECT_TRUE(observer->error);
  EXPECT_NE(NEM::MSH::oshGeoMesh::SafeDownCast(sb->GetOutputDataObject(0)),
            nullptr);
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
