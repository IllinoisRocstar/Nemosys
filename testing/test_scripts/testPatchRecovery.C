#include <patchRecovery.H>
#include <gtest.h>
#include <chrono>

const char* nodeMesh;
const char* recoveredMesh;
const char* hexmesh;

TEST(PatchRecoveryConstructor, ConstructWithArray)
{
  // load reference node mesh
  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(nodeMesh);
  const std::vector<int> arrayIDs = {0,1,2,3,4,5,7};
  int order = 1;
  std::unique_ptr<PatchRecovery> recoverObj
    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh.get(),order,arrayIDs));
  recoverObj->computeNodalError();
  mesh->write("errorTest.vtu");
}

//TEST(PatchRecoveryTensorProdBasis, RecoverNodalSol)
//{
//  std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(hexmesh);
//  const std::vector<int> arrayIDs = {0,1};
//  int order = 1;
//  std::unique_ptr<PatchRecovery> recoverObj
//    = std::unique_ptr<PatchRecovery> (new PatchRecovery(mesh.get(),order,arrayIDs));
//  Timer T;
//  T.start();
//  recoverObj->recoverNodalSolution(0);
//  T.stop();
//  std::cout << "time for regpoly (ms): " << T.elapsed() << std::endl;
//  mesh->write("polyApproxRecovery.vtu");
//  T.start();
//  recoverObj->recoverNodalSolution(1);
//  T.stop();
//  std::cout << "time for orthopoly(ms): " << T.elapsed() << std::endl;
//  mesh->write("orthoPolyApproxRecovery.vtu");
//}
 
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 4);
  nodeMesh = argv[1];
  recoveredMesh = argv[2];
  hexmesh = argv[3];
  return RUN_ALL_TESTS();
}
