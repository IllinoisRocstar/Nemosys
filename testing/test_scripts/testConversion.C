#include <meshBase.H>
#include <foamMesh.H>
#include <gtest.h>

const char* mshName;
const char* volName;
const char* refMshVTUName;
const char* refVolVTUName;
const char* legacyVTK1;
const char* legacyVTK2;
const char* legacyVTK1_ref;
const char* legacyVTK2_ref;
const char* stlvtp;
const char* pnttri;
const char* pnttri_ref;
const char* pntquad;
const char* pntquad_ref;
const char* pnthex;
const char* pnthex_ref;
const char* pntmix;
const char* pntmix_ref;
const char* packConv;
const char* packConv_ref;

TEST(Conversion, ConvertGmshToVTK)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(mshName);
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(refMshVTUName);
  EXPECT_EQ(0,diffMesh(mesh.get(),refMesh.get())); 
} 

TEST(Conversion, ConvertVolToVTK)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(volName);
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(refVolVTUName);
  EXPECT_EQ(0,diffMesh(mesh.get(),refMesh.get()));
}

TEST(Conversion, ConvertLegacyVTKToVTU)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(legacyVTK1);
  std::unique_ptr<meshBase> mesh_ref = meshBase::CreateUnique(legacyVTK1_ref);
  std::unique_ptr<meshBase> mesh1 = meshBase::CreateUnique(legacyVTK2);
  std::unique_ptr<meshBase> mesh1_ref = meshBase::CreateUnique(legacyVTK2_ref);
  std::cout << mesh->getNumberOfCells() << std::endl;
  std::cout << mesh_ref->getNumberOfCells() << std::endl;
  std::cout << mesh1->getNumberOfCells() << std::endl;
  std::cout << mesh1_ref->getNumberOfCells() << std::endl;
  EXPECT_EQ(0, diffMesh(mesh.get(), mesh_ref.get()));
  EXPECT_EQ(0, diffMesh(mesh1.get(), mesh1_ref.get()));
}

TEST(Conversion, ConvertVTPToSTLAndBack)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(stlvtp);
  mesh->write("gorilla-test.stl");
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique("gorilla-test.stl");
  EXPECT_EQ(0, diffMesh(mesh.get(), refMesh.get()));
  if (remove("gorilla-test.stl"))
  {
    std::cerr << "Error removing gorilla-test.stl" << std::endl;
    exit(1);
  } 
}

TEST(Conversion, ConvertPNTTriToVTU)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(pnttri);
  mesh->write("pnt-tri-test.vtu");
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(pnttri_ref);
  EXPECT_EQ(0, diffMesh(mesh.get(), refMesh.get()));
  if (remove("pnt-tri-test.vtu"))
  {
    std::cerr << "Error removing pnt-tri-test.vtu" << std::endl;
    exit(1);
  } 
}

TEST(Conversion, ConvertPNTQuadToVTU)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(pntquad);
  mesh->write("pnt-quad-test.vtu");
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(pntquad_ref);
  EXPECT_EQ(0, diffMesh(mesh.get(), refMesh.get()));
  if (remove("pnt-quad-test.vtu"))
  {
    std::cerr << "Error removing pnt-quad-test.vtu" << std::endl;
    exit(1);
  } 
}

TEST(Conversion, ConvertPNTHexToVTU)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(pnthex);
  mesh->write("pnt-hex-test.vtu");
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(pnthex_ref);
  EXPECT_EQ(0, diffMesh(mesh.get(), refMesh.get()));
  if (remove("pnt-hex-test.vtu"))
  {
    std::cerr << "Error removing pnt-hex-test.vtu" << std::endl;
    exit(1);
  } 
}

TEST(Conversion, ConvertPNTMixToVTU)
{
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(pntmix);
  mesh->write("pnt-mix-test.vtu");
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(pntmix_ref);
  EXPECT_EQ(0, diffMesh(mesh.get(), refMesh.get()));
  if (remove("pnt-mix-test.vtu"))
  {
    std::cerr << "Error removing pnt-mix-test.vtu" << std::endl;
    exit(1);
  } 
}

#ifdef HAVE_CFMSH
TEST(Conversion, ConvertVTUToFoam)
{
	// create meshBase object
  std::shared_ptr<meshBase> mesh = meshBase::CreateShared(packConv);

  // create foamMesh object
  FOAM::foamMesh *fm = new FOAM::foamMesh(mesh);

  // Write polyMesh
  fm->write("name");

  // Reads polyMesh
  meshBase *fm2 = new FOAM::foamMesh();
  fm2->read("name");

  std::unique_ptr<meshBase> origMesh = meshBase::CreateUnique(fm2);

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(packConv_ref);

  EXPECT_EQ(0, diffMesh(origMesh.get(), refMesh.get()));
}
#endif

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 20);
  refMshVTUName = argv[1];
  mshName = argv[2];
  refVolVTUName = argv[3];
  volName = argv[4];
  legacyVTK1 = argv[5];
  legacyVTK2 = argv[6];
  legacyVTK1_ref = argv[7];
  legacyVTK2_ref = argv[8];
  stlvtp = argv[9];
  pnttri = argv[10];
  pnttri_ref = argv[11];
  pntquad = argv[12];
  pntquad_ref = argv[13];
  pnthex = argv[14];
  pnthex_ref = argv[15];
  pntmix = argv[16];
  pntmix_ref = argv[17];
  packConv = argv[18];
  packConv_ref = argv[19];
  return RUN_ALL_TESTS();
}

