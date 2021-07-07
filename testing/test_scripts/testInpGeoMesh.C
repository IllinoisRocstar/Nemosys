#include <gtest/gtest.h>
#include <Mesh/diffMesh.H>
#include <Mesh/geoMeshFactory.H>
#include <Mesh/inpGeoMesh.H>

namespace {
std::string gmshMesh;
std::string inpMesh;
std::string vtkMesh;
}  // namespace

TEST(inpGeoMesh, gmsh2Inp2GM) {
  auto gmshGM =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(gmshMesh));
  vtkNew<NEM::MSH::inpGeoMesh> inpGM;
  inpGM->takeGeoMesh(gmshGM);
  auto inpFile =
      "test_" + gmshMesh.substr(0, gmshMesh.find_last_of('.')) + ".inp";
  inpGM->write(inpFile);
  auto reRead =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(inpFile));
  EXPECT_EQ(0, NEM::MSH::diffMesh(inpGM, reRead, 1e-6, 1e-4));
}

TEST(inpGeoMesh, inp2Vtk) {
  auto inpGM =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(inpMesh));
  auto vtkGM = vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(
      NEM::MSH::New(NEM::MSH::MeshType::VTK_GEO_MESH));
  vtkGM->takeGeoMesh(inpGM);
  vtkGM->write("test_" + inpMesh.substr(0, inpMesh.find_last_of('.')) + ".vtu");
  auto gold =
      vtkSmartPointer<NEM::MSH::geoMeshBase>::Take(NEM::MSH::Read(vtkMesh));
  EXPECT_EQ(0, NEM::MSH::diffMesh(vtkGM, gold));
}

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  assert(argc == 4);
  gmshMesh = argv[1];
  inpMesh = argv[2];
  vtkMesh = argv[3];

  int res = RUN_ALL_TESTS();

  return res;
}
