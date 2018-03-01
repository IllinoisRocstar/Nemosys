#include <patchRecovery.H>
#include <gtest.h>

const char* nodeMesh;


TEST(PatchRecoveryConstructor, ConstructWithoutArray)
{
  // load reference node mesh
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(nodeMesh);
  const std::vector<int> arrayIDs = {0,2,3,7};
  int order = 2;
  std::unique_ptr<PatchRecovery> recoverObj
    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh.get(),order,arrayIDs));
  recoverObj->recoverNodalSolution();
}

 
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 2);
  nodeMesh = argv[1];
  return RUN_ALL_TESTS();
}
