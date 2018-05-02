#include <OrderOfAccuracy.H>
#include <gtest.h>

const char* coarse;
const char* fine;
const char* finer;
const char* rich_ref;

// The fixture for testing class orthoPoly. From google test primer.
class OrderOfAccuracyTest : public ::testing::Test 
{
  protected:
    // You can remove any or all of the following functions if its body
    // is empty.
  
    OrderOfAccuracyTest()
    {
      f1 = meshBase::CreateUnique(finer);
      f2 = meshBase::CreateUnique(fine);
      f3 = meshBase::CreateUnique(coarse);
      oac = new OrderOfAccuracy(f1.get(),f2.get(),f3.get(),arrayIDs);
    }
  
    virtual ~OrderOfAccuracyTest()
    {}
  
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:
    virtual void SetUp()
    {
  
    }
  
    virtual void TearDown() {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }
  
    // Objects declared here can be used by all tests in the test case for orthoPoly.
    std::unique_ptr<meshBase> f1;
    std::unique_ptr<meshBase> f2;
    std::unique_ptr<meshBase> f3;
    const std::vector<int> arrayIDs = {0};
    const std::vector<std::vector<double>> res_ref = {{0.32581996578105021,
                                                       0.19338485350409632,
                                                       0.32030546671960047}};
    const std::vector<std::vector<double>> oac_ref = {{2.1119174483287679,
                                                       1.9704754800492528,
                                                       2.0774696352142827}};
    const std::vector<std::vector<double>> asymp_ref = {{1.0200502831373015,
                                                         1.0390518032820475,
                                                         1.031263791985966}};
    OrderOfAccuracy* oac; 
};

// tests richardson extrapolation 
TEST_F(OrderOfAccuracyTest, ComputeRichardsonExtrap)
{
  oac->computeRichardsonExtrapolation();
  std::unique_ptr<meshBase> f1_ref = meshBase::CreateUnique(rich_ref);
  EXPECT_EQ(0,diffMesh(f1.get(),f1_ref.get()));
}

// tests if correct resoltution desired error is returned
TEST_F(OrderOfAccuracyTest, ComputeResolution)
{
  std::vector<std::vector<double>> res = oac->computeResolution(.005);
  for (int j = 0; j < res[0].size(); ++j)
  {
    EXPECT_DOUBLE_EQ(res_ref[0][j],res[0][j]);
  }

} 

// tests if correct oac is returned
TEST_F(OrderOfAccuracyTest, ComputeOrderOfAccuracy)
{
  std::vector<std::vector<double>> order = oac->computeOrderOfAccuracy(); 
  for (int j = 0; j < order[0].size(); ++j)
  {
    EXPECT_DOUBLE_EQ(order[0][j],oac_ref[0][j]);
  }
}

// tests if correct asymptotic ratios are returned
TEST_F(OrderOfAccuracyTest, CheckAsymptoticRange)
{
  std::vector<std::vector<double>> asymp = oac->checkAsymptoticRange();
  for (int j = 0; j < asymp[0].size(); ++j)
  {
    EXPECT_DOUBLE_EQ(asymp[0][j],asymp_ref[0][j]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  finer = argv[1];
  fine = argv[2];
  coarse = argv[3];
  rich_ref = argv[4];
  return RUN_ALL_TESTS();
}
