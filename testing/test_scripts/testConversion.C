#include <meshBase.H>
#include <gtest.h>

const char* mshName;
const char* volName;
const char* refMshVTUName;
const char* refVolVTUName;
const char* legacyVTK1;
const char* legacyVTK2;
const char* legacyVTK1_ref;
const char* legacyVTK2_ref;

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

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 9);
  refMshVTUName = argv[1];
  mshName = argv[2];
  refVolVTUName = argv[3];
  volName = argv[4];
  legacyVTK1 = argv[5];
  legacyVTK2 = argv[6];
  legacyVTK1_ref = argv[7];
  legacyVTK2_ref = argv[8];
  return RUN_ALL_TESTS();
}

