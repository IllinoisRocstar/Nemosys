#include <patchRecovery.H>
#include <gtest.h>

const char* nodeMesh;


TEST(PatchRecoveryConstructor, ConstructWithoutArray)
{
  // load reference node mesh
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(nodeMesh);
  const std::vector<int> arrayIDs = {0,1,2,3,4,5};
  int order = 1;
  std::unique_ptr<PatchRecovery> recoverObj
    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh.get(),order,arrayIDs));
  recoverObj->recoverNodalSolution();
  mesh->write("test.vtu");
}

 
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 2);
  nodeMesh = argv[1];
  return RUN_ALL_TESTS();
}
