#include <Mesh/geoMeshBase.H>
#include <Mesh/geoMeshFactory.H>

#include <Mesh/diffMesh.H>
#include <gtest/gtest.h>

std::string foam_name;
std::string vtk_ref_foam;
std::string vtk_name;
std::string foam_ref_vtk;
std::string gmsh_name;
std::string vtk_ref_gmsh;
std::string gmsh_ref;
std::string exo_name;
std::string vtk_ref_exo;
std::string exo_ref;
std::string osh_name;
std::string vtk_ref_osh;
std::string osh_ref;
std::string inp_name;
std::string vtk_ref_inp;
std::string inp_ref;

#ifdef HAVE_CFMSH
// convert foam to vtk
TEST(SmartConversion, ConvertFOAMToVTK) {
  // read foam mesh as foam
  NEM::MSH::geoMeshBase *foam_mesh = NEM::MSH::Read(foam_name, "foam");

  // read vtk reference mesh
  NEM::MSH::geoMeshBase *vtk_ref_mesh = NEM::MSH::Read(vtk_ref_foam);

  // new vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::New("vtk_mesh", "vtk");

  // convert
  vtk_mesh->takeGeoMesh(foam_mesh);

  EXPECT_EQ(0, NEM::MSH::diffMesh(vtk_mesh, vtk_ref_mesh));
}

// convert vtk to foam
TEST(SmartConversion, ConvertVTKToFOAM) {
  // read vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::Read(vtk_name, "vtk");

  // read foam reference mesh
  NEM::MSH::geoMeshBase *foam_ref_mesh = NEM::MSH::Read(foam_ref_vtk, "foam");

  // new foam mesh as foam
  NEM::MSH::geoMeshBase *foam_mesh = NEM::MSH::New("foam_mesh", "foam");

  // convert
  foam_mesh->takeGeoMesh(vtk_mesh);

  // foam meshes are changed when they are written--this may fail if only
  // compared in memory.
  foam_mesh->write("test_foam_mesh");
  NEM::MSH::geoMeshBase *written_foam_mesh =
      NEM::MSH::Read("test_foam_mesh", "foam");

  EXPECT_EQ(0, NEM::MSH::diffMesh(written_foam_mesh, foam_ref_mesh));
}
#endif

#ifdef HAVE_GMSH
// convert gmsh to vtk
TEST(SmartConversion, ConvertGMSHToVTK) {
  // read gmsh mesh as gmsh
  NEM::MSH::geoMeshBase *gmsh_mesh = NEM::MSH::Read(gmsh_name, "gmsh");

  // read vtk reference mesh
  NEM::MSH::geoMeshBase *vtk_ref_mesh = NEM::MSH::Read(vtk_ref_gmsh, "vtk");

  // new vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::New("vtk_mesh", "vtk");

  // convert
  vtk_mesh->takeGeoMesh(gmsh_mesh);

  EXPECT_EQ(0, NEM::MSH::diffMesh(vtk_mesh, vtk_ref_mesh));
}

// convert vtk to gmsh
TEST(SmartConversion, ConvertVTKToGMSH) {
  // read vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::Read(vtk_name, "vtk");

  // read gmsh reference mesh
  NEM::MSH::geoMeshBase *gmsh_ref_mesh = NEM::MSH::Read(gmsh_ref);

  // new gmsh mesh as gmsh
  NEM::MSH::geoMeshBase *gmsh_mesh = NEM::MSH::New("gmsh_mesh", "gmsh");

  // convert
  gmsh_mesh->takeGeoMesh(vtk_mesh);

  EXPECT_EQ(0, NEM::MSH::diffMesh(gmsh_mesh, gmsh_ref_mesh));
}
#endif

// convert exodus ii to vtk
TEST(SmartConversion, ConvertEXOToVTK) {
  // read exo mesh as exo
  NEM::MSH::geoMeshBase *exo_mesh = NEM::MSH::Read(exo_name, "exodus");

  // read vtk reference mesh
  NEM::MSH::geoMeshBase *vtk_ref_mesh = NEM::MSH::Read(vtk_ref_exo);

  // new vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::New("vtk_mesh", "vtk");

  // convert
  vtk_mesh->takeGeoMesh(exo_mesh);

  EXPECT_EQ(0, NEM::MSH::diffMesh(vtk_mesh, vtk_ref_mesh));
}

// convert vtk to exodus ii
TEST(SmartConversion, ConvertVTKToEXO) {
  // read vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::Read(vtk_name, "vtk");

  // read exodus reference mesh
  NEM::MSH::geoMeshBase *exo_ref_mesh = NEM::MSH::Read(exo_ref);

  // new exodus mesh as exodus
  NEM::MSH::geoMeshBase *exo_mesh = NEM::MSH::New("exo_mesh", "exodus");

  // convert
  exo_mesh->takeGeoMesh(vtk_mesh);

  // exodus meshes are changed when they are written--this may fail if only
  // compared in memory.
  exo_mesh->write("test_exo_mesh.exo");
  NEM::MSH::geoMeshBase *written_exo_mesh =
      NEM::MSH::Read("test_exo_mesh.exo", "exodus");

  EXPECT_EQ(0, NEM::MSH::diffMesh(written_exo_mesh, exo_ref_mesh));
}

// convert omega h to vtk
TEST(SmartConversion, ConvertOSHToVTK) {
  // read osh mesh as osh
  NEM::MSH::geoMeshBase *osh_mesh = NEM::MSH::Read(osh_name, "omegah");

  // read vtk reference mesh
  NEM::MSH::geoMeshBase *vtk_ref_mesh = NEM::MSH::Read(vtk_ref_osh);

  // new vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::New("vtk_mesh", "vtk");

  // convert
  vtk_mesh->takeGeoMesh(osh_mesh);

  EXPECT_EQ(0, NEM::MSH::diffMesh(vtk_mesh, vtk_ref_mesh));
}

// convert vtk to omega h
TEST(SmartConversion, ConvertVTKToOSH) {
  // read vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::Read(vtk_name, "vtk");

  // read omegah reference mesh
  NEM::MSH::geoMeshBase *osh_ref_mesh = NEM::MSH::Read(osh_ref);

  // new omegah mesh as omegah
  NEM::MSH::geoMeshBase *osh_mesh = NEM::MSH::New("osh_mesh", "omegah");

  // convert
  osh_mesh->takeGeoMesh(vtk_mesh);

  EXPECT_EQ(0, NEM::MSH::diffMesh(osh_mesh, osh_ref_mesh));
}

// convert calculix/abaqus to vtk
TEST(SmartConversion, ConvertINPToVTK) {
  // read inp mesh as inp
  NEM::MSH::geoMeshBase *inp_mesh = NEM::MSH::Read(inp_name, "calculix");

  // read vtk reference mesh
  NEM::MSH::geoMeshBase *vtk_ref_mesh = NEM::MSH::Read(vtk_ref_inp);

  // new vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::New("vtk_mesh", "vtk");

  // convert
  vtk_mesh->takeGeoMesh(inp_mesh);

  EXPECT_EQ(0, NEM::MSH::diffMesh(vtk_mesh, vtk_ref_mesh));
}

// convert vtk to calculix/abaqus
TEST(SmartConversion, ConvertVTKToINP) {
  // read vtk mesh as vtk
  NEM::MSH::geoMeshBase *vtk_mesh = NEM::MSH::Read(vtk_name, "vtk");

  // read inp reference mesh
  NEM::MSH::geoMeshBase *inp_ref_mesh = NEM::MSH::Read(inp_ref);

  // new inp mesh as inp
  NEM::MSH::geoMeshBase *inp_mesh = NEM::MSH::New("inp_mesh", "calculix");

  // convert
  inp_mesh->takeGeoMesh(vtk_mesh);

  // inp meshes are changed when they are written--this may fail if only
  // compared in memory.
  inp_mesh->write("test_inp_mesh.inp");
  NEM::MSH::geoMeshBase *written_inp_mesh =
      NEM::MSH::Read("test_inp_mesh.inp", "calculix");

  EXPECT_EQ(0, NEM::MSH::diffMesh(written_inp_mesh, inp_ref_mesh));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 17);
  foam_name = argv[1];
  vtk_ref_foam = argv[2];
  vtk_name = argv[3];
  foam_ref_vtk = argv[4];
  gmsh_name = argv[5];
  vtk_ref_gmsh = argv[6];
  gmsh_ref = argv[7];
  exo_name = argv[8];
  vtk_ref_exo = argv[9];
  exo_ref = argv[10];
  osh_name = argv[11];
  vtk_ref_osh = argv[12];
  osh_ref = argv[13];
  inp_name = argv[14];
  vtk_ref_inp = argv[15];
  inp_ref = argv[16];
  return RUN_ALL_TESTS();
}
