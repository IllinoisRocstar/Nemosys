#include <Cubature.H>
#include <gtest.h>

const char* nodeMesh;
const char* refGauss;
// The fixture for testing class orthoPoly. From google test primer.
class CubatureTest : public ::testing::Test 
{
  protected:
    // You can remove any or all of the following functions if its body
    // is empty.
  
    CubatureTest()
    {
      mesh = meshBase::Create(nodeMesh);
      cuby = new GaussCubature(mesh);
		}
  
    virtual ~CubatureTest()
    {
      if (mesh) 
      {
        delete mesh;
        mesh = 0;
      }
      if (cuby) 
      {
        delete cuby;
        cuby = 0;
      }
    }
  
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
    meshBase* mesh;
    GaussCubature* cuby; 
    const std::vector<int> arrayIDs = {0,2,3,7}; 
};


// tests construction, writeout, loading and interpolation
TEST_F(CubatureTest, InterpolateToGauss)
{
  cuby->constructGaussMesh(arrayIDs);
  meshBase* gaussMesh = meshBase::Create("gaussTest.vtp");
  meshBase* refGaussMesh = meshBase::Create(refGauss);
  EXPECT_EQ(0,diffMesh(gaussMesh,refGaussMesh));
	if (gaussMesh) delete gaussMesh;
	if (refGaussMesh) delete refGaussMesh;
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 3);
  nodeMesh = argv[1];
  refGauss = argv[2];
  return RUN_ALL_TESTS();
}

