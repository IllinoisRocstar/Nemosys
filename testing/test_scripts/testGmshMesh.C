#include <gmshMesh.H>
#include <gtest.h>

#include <sstream>

std::string arg_fnameMSH2;
std::string arg_fnameMSH4;
std::string arg_fnameMSH41;

// Aux functions
bool filesEqual(const std::string &fn1, const std::string &fn2)
{
  std::fstream f1, f2;
  f1.open(fn1, fstream::in);
  f2.open(fn2, fstream::in);
  char string1[256], string2[256];
  int j = 0;
  while (!f1.eof()) {
    f1.getline(string1, 256);
    f2.getline(string2, 256);
    j++;
    if (j == 2) {
      std::stringstream ss1(string1);
      std::stringstream ss2(string2);
      float ver1, ver2;
      int bin1, bin2, bit1, bit2;
      ss1 >> ver1 >> bin1 >> bit1;
      ss2 >> ver2 >> bin2 >> bit2;
      if ((bit1 == 8 && bit2 == 4)
          || (bit1 == 4 && bit2 == 8))
      {
        std::cout << "ERROR: Underlying Gmsh is built with floats."
                  << " NEMoSys requires support for doubles."
                  << "\n";
          return false;
      }
    }
    if (strcmp(string1, string2) != 0) {
      std::cout << j << "-the strings are not equal" << "\n";
      std::cout << "   " << string1 << "\n";
      std::cout << "   " << string2 << "\n";
      return false;
    }
  }
  return true;
}

int test_gmshMeshDefaultConstructor() {
  gmshMesh gmshMesh1 = gmshMesh();
  return 0;
}

int test_gmshMeshFilenameConstructor(const std::string &fname) {
  auto *gmshMesh1 = new gmshMesh(fname);

  if (gmshMesh1->getFileName() != fname) return 1;

  if (gmshMesh1->getNumberOfPoints() != 177) return 1;
  if (gmshMesh1->getNumberOfCells() != 573) return 1;

  delete gmshMesh1;

  return 0;
}

int test_gmshMeshWriteMSH2(const std::string &fname) {
  std::string test_fname = fname + "_test";
  gmshMesh gmshMesh1 = gmshMesh(fname);

  gmshMesh1.write(test_fname, 2.0, false);

  return filesEqual(fname, test_fname);
}

int test_gmshMeshWriteMSH4(const std::string &fname) {
  std::string test_fname = fname + "_test";
  gmshMesh gmshMesh1 = gmshMesh(fname);

  gmshMesh1.write(test_fname, 4.0, false);

  return filesEqual(fname, test_fname);
}

int test_gmshMeshWriteMSH41(const std::string &fname) {
  std::string test_fname = fname + "_test";
  gmshMesh gmshMesh1 = gmshMesh(fname);

  gmshMesh1.write(test_fname, 4.1, false);

  return filesEqual(fname, test_fname);
}

// TEST macros
TEST(gmshMesh, DefaultConstructor)
{
  EXPECT_EQ(0, test_gmshMeshDefaultConstructor());
}

TEST(gmshMesh, FilenameConstructor)
{
  EXPECT_EQ(0, test_gmshMeshFilenameConstructor(arg_fnameMSH41));
}

TEST(gmshMesh, WriteMSH2)
{
  EXPECT_EQ(1, test_gmshMeshWriteMSH2(arg_fnameMSH2));
}

TEST(gmshMesh, WriteMSH4)
{
  EXPECT_EQ(1, test_gmshMeshWriteMSH4(arg_fnameMSH4));
}

TEST(gmshMesh, WriteMSH41)
{
  EXPECT_EQ(1, test_gmshMeshWriteMSH41(arg_fnameMSH41));
}

// test constructor
int main(int argc, char** argv) {
  // IO
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 4);
  arg_fnameMSH2 = argv[1];
  arg_fnameMSH4 = argv[2];
  arg_fnameMSH41 = argv[3];

  // run tests
  int res = RUN_ALL_TESTS();

  return res;
}
