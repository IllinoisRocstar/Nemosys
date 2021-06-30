#include <exoMesh.H>
#include <foamMesh.H>
#include <meshBase.H>
#include <Drivers/NemDriver.H>
#include <Drivers/Conversion/ManipExoConversionDriver.H>
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
const char* buildingTet;
const char* buildingTet_ref;
const char* combineBlocks_json;

std::string combine_blocks(const char *jsonF) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }
  jsoncons::json inputjson;
  inputStream >> inputjson;

  auto nemdrvobj = NEM::DRV::NemDriver::readJSON(inputjson);
  EXPECT_NE(nemdrvobj, nullptr);
  nemdrvobj->execute();

  std::string ofname =
      dynamic_cast<NEM::DRV::ManipExoConversionDriver *>(nemdrvobj.get())
          ->getFiles()
          .outputMeshFile;
  return ofname;
}

TEST(Conversion, ConvertEXOtoEXO) {
  auto mesh = combine_blocks(combineBlocks_json);
  auto em = new NEM::MSH::EXOMesh::exoMesh(mesh);
  em->read(mesh);

  int numBlks = em->getNumberOfElementBlocks();
  std::cout << "Number of element blocks " << numBlks << std::endl;
  int numSdes = em->getNumSdesInSdeSetById(4);
  std::cout << "Number of sides in sideset " << numSdes << std::endl;

  int ret;
  if (numBlks != 3 || numSdes != 52) {
    std::cout << "Expected number of blocks to be 3 and sides to be 52."
              << std::endl;
    ret = 1;
  } else
    ret = 0;
  EXPECT_EQ(0, ret);
}

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

TEST(Conversion, ConvertVTKHexToTet)
{
  // create meshBase object
  std::shared_ptr<meshBase> mesh = meshBase::CreateShared(buildingTet);
  
  // Converts hex mesh to tet mesh and writes in VTU file.
  mesh->convertHexToTetVTK(mesh->getDataSet());
  mesh->write("outputTet.vtu");

  std::unique_ptr<meshBase> origMesh = meshBase::CreateUnique("outputTet.vtu");

  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(buildingTet_ref);

  EXPECT_EQ(0, diffMesh(origMesh.get(), refMesh.get()));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 23);
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
  buildingTet = argv[20];
  buildingTet_ref = argv[21];
  combineBlocks_json = argv[22];
  return RUN_ALL_TESTS();
}

