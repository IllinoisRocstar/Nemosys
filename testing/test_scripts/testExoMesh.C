#include <gtest/gtest.h>
#include <Mesh/exoMesh.H>

std::string arg_fName1;
std::string arg_fName2;
std::string arg_fName3;
std::string arg_fName4;
std::string four_boxes_tetra;
std::string four_boxes_hexa;
std::string simple_circles;

bool diffExoMesh(const NEM::MSH::EXOMesh::exoMesh &em1,
                 const NEM::MSH::EXOMesh::exoMesh &em2) {
  if (em1.getNumberOfNodes() != em2.getNumberOfNodes()) return true;
  if (em1.getNumberOfElements() != em2.getNumberOfElements()) return true;
  if (em1.getNumberOfNodeSets() != em2.getNumberOfNodeSets()) return true;
  if (em1.getNumberOfElementBlocks() != em2.getNumberOfElementBlocks())
    return true;
  if (em1.getNumberOfSideSets() != em2.getNumberOfSideSets()) return true;

  return false;
}

int test_exoMesh_stitch(const std::string &fName1, const std::string &fName2) {
  NEM::MSH::EXOMesh::exoMesh em1;
  NEM::MSH::EXOMesh::exoMesh em2;

  em1.read(fName1);
  em2.read(fName2);

  int non = em1.getNumberOfNodes() + em2.getNumberOfNodes();
  int noe = em1.getNumberOfElements() + em2.getNumberOfElements();
  int nons = em1.getNumberOfNodeSets() + em2.getNumberOfNodeSets();
  int noeb = em1.getNumberOfElementBlocks() + em2.getNumberOfElementBlocks();
  int noss = em1.getNumberOfSideSets() + em2.getNumberOfSideSets();

  em1.stitch(em2);

  em1.setFileName("test_stitch.g");
  em1.write();
  em1.report();

  if (em1.getNumberOfNodes() != non) return 1;
  if (em1.getNumberOfElements() != noe) return 1;
  if (em1.getNumberOfNodeSets() != nons) return 1;
  if (em1.getNumberOfElementBlocks() != noeb) return 1;
  if (em1.getNumberOfSideSets() != noss) return 1;

  return 0;
}

bool test_exoMesh_mergeNodes(const std::string &fName1,
                             const std::string &fName2) {
  NEM::MSH::EXOMesh::exoMesh em1;
  NEM::MSH::EXOMesh::exoMesh em2;

  em1.read(fName1);
  em2.read(fName2);

  em1.mergeNodes();

  em1.setFileName("test_merge.g");
  em1.write();
  em1.report();

  return !diffExoMesh(em1, em2);
}

int test_exoMesh_combineBlocks_tetra(const std::string &fname) {
  auto em = new NEM::MSH::EXOMesh::exoMesh(fname);
  em->read(fname);

  em->combineElmBlks({1, 2}, "");

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

  return ret;
}

int test_exoMesh_combineBlocks_hexa(const std::string &fname) {
  auto em = new NEM::MSH::EXOMesh::exoMesh(fname);
  em->read(fname);

  em->combineElmBlks({1, 2}, "");
  em->combineElmBlks({3, 4}, "");

  int numBlks = em->getNumberOfElementBlocks();
  std::cout << "Number of element blocks " << numBlks << std::endl;
  int numSdes = em->getNumSdesInSdeSetById(3);
  std::cout << "Number of sides in sideset " << numSdes << std::endl;

  int ret;
  if (numBlks != 2 || numSdes != 50) {
    std::cout << "Expected number of blocks to be 2 and sides to be 50."
              << std::endl;
    ret = 1;
  } else
    ret = 0;

  return ret;
}

int test_exoMesh_combineBlocks_circles(const std::string &fname) {
  auto em = new NEM::MSH::EXOMesh::exoMesh(fname);
  em->read(fname);

  em->combineElmBlks({1, 2}, "");

  int numBlks = em->getNumberOfElementBlocks();
  std::cout << "Number of element blocks " << numBlks << std::endl;
  int numSdes = em->getNumSdesInSdeSetById(1);
  std::cout << "Number of sides in sideset " << numSdes << std::endl;

  int ret;
  if (numBlks != 1 || numSdes != 52) {
    std::cout << "Expected number of blocks to be 1 and sides to be 52."
              << std::endl;
    ret = 1;
  } else
    ret = 0;

  return ret;
}

TEST(exoMesh, stitch) {
  EXPECT_EQ(0, test_exoMesh_stitch(arg_fName1, arg_fName2));
}

TEST(exoMesh, mergeNodes) {
  EXPECT_EQ(true, test_exoMesh_mergeNodes(arg_fName3, arg_fName4));
}

TEST(exoMesh, combineBlocks_tetra) {
  EXPECT_EQ(0, test_exoMesh_combineBlocks_tetra(four_boxes_tetra));
}

TEST(exoMesh, combineBlocks_hexa) {
  EXPECT_EQ(0, test_exoMesh_combineBlocks_hexa(four_boxes_hexa));
}

TEST(exoMesh, combineBlocks_circles) {
  EXPECT_EQ(0, test_exoMesh_combineBlocks_circles(simple_circles));
}

// test constructor
int main(int argc, char **argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 8);
  arg_fName1 = argv[1];
  arg_fName2 = argv[2];
  arg_fName3 = argv[3];
  arg_fName4 = argv[4];
  four_boxes_tetra = argv[5];
  four_boxes_hexa = argv[6];
  simple_circles = argv[7];

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
