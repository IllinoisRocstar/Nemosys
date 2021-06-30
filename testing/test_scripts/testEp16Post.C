#include <gtest/gtest.h>
#include <cmath>
#include "ep16Post.H"
#include "point.H"

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
  std::vector<Point> means(ep->kmeans->getMeans());
  ASSERT_EQ(4, means.size());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
