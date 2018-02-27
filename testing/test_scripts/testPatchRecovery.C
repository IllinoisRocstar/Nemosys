#include <patchRecovery.H>
#include <gtest.h>

const char* nodeMesh;
const char* recoveredMesh;

TEST(PatchRecoveryConstructor, ConstructWithoutArray)
{
  // load reference node mesh
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(nodeMesh);
  const std::vector<int> arrayIDs = {0,1,2,3,4,5,7};
  int order = 1;
  std::unique_ptr<PatchRecovery> recoverObj
    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh.get(),order,arrayIDs));
  recoverObj->recoverNodalSolution();
  std::unique_ptr<meshBase> refMesh = meshBase::CreateUnique(recoveredMesh);
  EXPECT_EQ(0,diffMesh(mesh.get(),refMesh.get()));
}

 
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 3);
  nodeMesh = argv[1];
  recoveredMesh = argv[2];
  return RUN_ALL_TESTS();
}
