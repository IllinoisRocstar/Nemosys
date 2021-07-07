#include <Drivers/Conversion/ConversionDriver.H>
#include <gtest/gtest.h>
const char *vtk2patranJSON;
const char *vtk2patran_out;
const char *vtk2patran_ref;

// Test VTK->PATRAN conversion driver
// Note quite crude - tests for equality by character in file
// except file name (2nd line) and date (4th line)
TEST(Conversion, ConvertVTKToPatran) {
  std::string fName(vtk2patranJSON);
  std::ifstream inputStream(fName);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << vtk2patranJSON << std::endl;
  }
  jsoncons::json inputjson;
  inputStream >> inputjson;

  auto driver = NEM::DRV::NemDriver::readJSON(inputjson[0]);
  EXPECT_NE(driver, nullptr);
  driver->execute();
  std::ifstream resultFile(vtk2patran_out);
  EXPECT_TRUE(resultFile.good());
  std::ifstream reference(vtk2patran_ref);
  EXPECT_TRUE(reference.good());
  {
    std::string lineRef, lineResult;
    // 2nd line is file name and 4th line is date, so ignore them
    for (int i = 0; i < 4; ++i) {
      EXPECT_TRUE(std::getline(reference, lineRef));
      EXPECT_TRUE(std::getline(resultFile, lineResult));
      if (i % 2 == 0) {
        EXPECT_EQ(lineRef, lineResult);
      }
    }
    while (std::getline(reference, lineRef)) {
      EXPECT_TRUE(std::getline(resultFile, lineResult));
      EXPECT_EQ(lineRef, lineResult);
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 4);
  vtk2patranJSON = argv[1];
  vtk2patran_out = argv[2];
  vtk2patran_ref = argv[3];
  return RUN_ALL_TESTS();
}
