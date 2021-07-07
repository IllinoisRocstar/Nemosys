#include <gtest/gtest.h>

#include <InputGeneration/ep16Post.H>
#include <vtkSTLReader.h>
#include <string>

// data loading
class TestDriver : public testing::Test {
 protected:
  jsoncons::json inp;

  virtual void SetUp() {
    std::ifstream inputStream;
    inputStream.open("ep16post.json");
    inputStream >> inp;
  }
};

TEST_F(TestDriver, EndToEnd) {
  int ret = 0;
  NEM::EPC::ep16Post *ep =
      NEM::EPC::ep16Post::readJSON(inp.at("Input Generation Options"), ret);
  ASSERT_EQ(0, ret);
  delete ep;
  vtkNew<vtkSTLReader> stlReader;
  for (int i = 0; i < 4; ++i) {
    std::string stlFileName{"clust_" + std::to_string(i) + ".stl"};
    stlReader->SetFileName(stlFileName.c_str());
    stlReader->Update();
    ASSERT_NE(stlReader->GetOutput(), nullptr);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
