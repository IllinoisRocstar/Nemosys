#include <gtest.h>

#include <NemDriver.H>
#include <exoMesh.H>
#include <meshBase.H>

std::string arg_fName1;
std::string arg_fName2;
std::string arg_fName3;
std::string arg_fName4;

bool diffExoMesh(const EXOMesh::exoMesh &em1, const EXOMesh::exoMesh &em2) {
  if (em1.getNumberOfNode() != em2.getNumberOfNode()) return true;
  if (em1.getNumberOfElement() != em2.getNumberOfElement()) return true;
  if (em1.getNumberOfNodeSet() != em2.getNumberOfNodeSet()) return true;
  if (em1.getNumberOfElementBlock() != em2.getNumberOfElementBlock())
    return true;
  if (em1.getNumberOfSideSets() != em2.getNumberOfSideSets()) return true;

  return false;
}

int test_exoMesh_stitch(const std::string &fName1, const std::string &fName2) {
  EXOMesh::exoMesh em1;
  EXOMesh::exoMesh em2;

  em1.read(fName1);
  em2.read(fName2);

  int non = em1.getNumberOfNode() + em2.getNumberOfNode();
  int noe = em1.getNumberOfElement() + em2.getNumberOfElement();
  int nons = em1.getNumberOfNodeSet() + em2.getNumberOfNodeSet();
  int noeb = em1.getNumberOfElementBlock() + em2.getNumberOfElementBlock();
  int noss = em1.getNumberOfSideSets() + em2.getNumberOfSideSets();

  em1.stitch(em2);

  em1.setFileName("test_stitch.g");
  em1.write();
  em1.report();

  if (em1.getNumberOfNode() != non) return 1;
  if (em1.getNumberOfElement() != noe) return 1;
  if (em1.getNumberOfNodeSet() != nons) return 1;
  if (em1.getNumberOfElementBlock() != noeb) return 1;
  if (em1.getNumberOfSideSets() != noss) return 1;

  return 0;
}

bool test_exoMesh_mergeNodes(const std::string &fName1,
                             const std::string &fName2) {
  EXOMesh::exoMesh em1;
  EXOMesh::exoMesh em2;

  em1.read(fName1);
  em2.read(fName2);

  em1.mergeNodes();

  em1.setFileName("test_merge.g");
  em1.write();
  em1.report();

  return !diffExoMesh(em1, em2);
}

TEST(exoMesh, stitch) {
  EXPECT_EQ(0, test_exoMesh_stitch(arg_fName1, arg_fName2));
}

TEST(exoMesh, mergeNodes) {
  EXPECT_EQ(true, test_exoMesh_mergeNodes(arg_fName3, arg_fName4));
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  arg_fName1 = argv[1];
  arg_fName2 = argv[2];
  arg_fName3 = argv[3];
  arg_fName4 = argv[4];

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
