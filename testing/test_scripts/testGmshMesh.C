#include <gmshMesh.H>
#include <gtest.h>

std::string arg_fname;

// Aux functions
int compareFiles(const std::string &fn1, const std::string &fn2)
{
  fstream f1, f2;
  int pass = 1;
  f1.open(fn1, fstream::in);
  f2.open(fn2, fstream::in);
  char string1[256], string2[256];
  int j = 0;
  while (!f1.eof()) {
    f1.getline(string1, 256);
    f2.getline(string2, 256);
    j++;
    if (strcmp(string1, string2) != 0) {
      cout << j << "-the strings are not equal" << "\n";
      cout << "   " << string1 << "\n";
      cout << "   " << string2 << "\n";
      pass = 0;
      break;
    }
  }
  return pass;
}

int test_gmshMeshDefaultConstructor() {
  gmshMesh gmshMesh1 = gmshMesh();
  return 0;
}

int test_gmshMeshFilenameConstructor(const std::string &fname) {
  auto *gmshMesh1 = new gmshMesh(fname);

  if (gmshMesh1->getFileName() != fname) return 1;
  if (gmshMesh1->getSFBool() != false) return 1;
  if (gmshMesh1->getOrder() != 1) return 1;

//  if (gmshMesh1.getNumberOfCells() != numCells) return 1;
//  if (gmshMesh1.getNumberOfPoints() != numPoints) return 1;
//  if (gmshMesh1.getDataSet())

  delete gmshMesh1;

  return 0;
}

int test_gmshMeshWrite(const std::string &fname) {
  std::string test_fname = "test_" + fname;
  gmshMesh gmshMesh1 = gmshMesh(fname);

  gmshMesh1.write(test_fname);

  return 0;
}

// TEST macros
TEST(gmshMesh, DefaultConstructor)
{
  EXPECT_EQ(0, test_gmshMeshDefaultConstructor());
}

TEST(gmshMesh, FilenameConstructor)
{
  EXPECT_EQ(0, test_gmshMeshFilenameConstructor(arg_fname));
}

TEST(gmshMesh, Write)
{
  EXPECT_EQ(0, test_gmshMeshWrite(arg_fname));
}

// test constructor
int main(int argc, char** argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 2);
  arg_fname = argv[1];

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
