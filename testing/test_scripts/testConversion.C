#include <meshBase.H>
#include <gtest.h>

const char* mshName;
const char* volName;
const char* refMshVTUName;
const char* refVolVTUName;

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


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  refMshVTUName = argv[1];
  mshName = argv[2];
  refVolVTUName = argv[3];
  volName = argv[4];
  return RUN_ALL_TESTS();
}

