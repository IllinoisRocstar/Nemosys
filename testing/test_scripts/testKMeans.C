#include <gtest/gtest.h>
#include "point.H"
#include "kmeans.H"
#include <cmath>

using namespace NEM::GEO;
using namespace NEM::MTH;

class TestPoint : public testing::Test 
{
  protected:
  Point p1;
  Point p2;

  virtual void SetUp() 
  {
    p1 = Point(1.0, 2.0, 3.0);
    p2 = Point(2.0, 3.0, 4.0);
  }
};

TEST_F(TestPoint, TestConstructorArguments) 
{
  Point p(4, 5, 6);
  EXPECT_EQ(3, p.dimensions_);
  EXPECT_EQ(-1, p.cluster_);
}

TEST_F(TestPoint, TestConstructorVector) 
{
  std::vector<double> vec{0, 1, 2, 3, 4, 5};
  Point point(vec);

  const int ref_size = (int) vec.size();

  EXPECT_EQ(ref_size, point.dimensions_);
  ASSERT_EQ(ref_size, point.data_.size());
  for (int idx = 0; idx < 6; ++idx)
    EXPECT_EQ((double) idx, point.data_[idx]);
}

TEST_F(TestPoint, TestUpdateCluster) 
{
  bool ret = p1.update(2);
  EXPECT_EQ(2, p1.cluster_);
  EXPECT_EQ(true, ret);
}

TEST_F(TestPoint, TestUpdateClusterNoChange) 
{
  bool ret = p1.update(2);
  EXPECT_TRUE(ret);
  ret = p1.update(2);
  EXPECT_EQ(2, p1.cluster_);
  EXPECT_FALSE(ret);
}

TEST_F(TestPoint, TestAdd) 
{
  p1.add(p2);
  EXPECT_EQ(3.0, p1.data_[0]);
  EXPECT_EQ(5.0, p1.data_[1]);
  EXPECT_EQ(7.0, p1.data_[2]);
}

TEST_F(TestPoint, TestDistance) 
{
  const double dist = Point::distance(p1, p2);
  EXPECT_EQ(sqrt(3), dist);
}


// data loading
class TestDataLoading : public testing::Test 
{
  protected:
  virtual void SetUp() {}
};

TEST_F(TestDataLoading, TestLoadA) 
{
  std::vector<Point> points;
  const std::string filepath = "test_data_3.txt";
  KMeans::loadPoints(filepath, &points);
  ASSERT_EQ(2, points.size());
  ASSERT_EQ(3, points[0].dimensions_);
  ASSERT_EQ(3, points[1].dimensions_);
  EXPECT_EQ(1.23, points[0].data_[0]);
  EXPECT_EQ(4.56, points[0].data_[1]);
  EXPECT_EQ(7.89, points[0].data_[2]);
}


// end-to-end
class TestEndToEnd : public testing::Test 
{
  protected:
  std::unique_ptr<KMeans> kmeans;

  virtual void SetUp() 
  {
    std::vector<Point> points;
    KMeans::loadPoints("data_200_2.txt", &points);

    kmeans.reset(new KMeans(2, 10));
    kmeans->init(points);
  }
};

TEST_F(TestEndToEnd, TestLoading) 
{
  ASSERT_EQ(200, kmeans->getPoints().size());
  ASSERT_EQ(2, kmeans->getPoints()[0].dimensions_);

  const std::vector<Point> means = kmeans->getMeans();
  ASSERT_EQ(2, means.size());

  // All means should be initialized to a random (non-zero) point in the data.
  for (const auto &mean : means) 
  {
    ASSERT_EQ(2, mean.dimensions_);
    for (double d : mean.data_)
      EXPECT_NE(0.0, d);
  }
}

TEST_F(TestEndToEnd, TestCompute)
{
  // Just check that we converge.
  bool ret = kmeans->run();
  ASSERT_TRUE(ret);
}


int main(int argc, char **argv) 
{
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

