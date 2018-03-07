#include <Cubature.H>
#include <gtest.h>

const char* nodeMesh;
const char* refGauss;
const char* refGauss1;
const char* refInt;

// The fixture for testing class orthoPoly. From google test primer.
class CubatureTest : public ::testing::Test 
{
  protected:
    // You can remove any or all of the following functions if its body
    // is empty.
  
    CubatureTest()
    {
      //mesh = meshBase::Create(nodeMesh);
      //cuby = new GaussCubature(mesh,arrayIDs);
			mesh = meshBase::CreateShared(nodeMesh);
      cuby = GaussCubature::CreateShared(mesh.get(), arrayIDs);
    }
  
    virtual ~CubatureTest()
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
    std::shared_ptr<meshBase> mesh;
    std::shared_ptr<GaussCubature> cuby; 
    const std::vector<int> arrayIDs = {0,2,3,7}; 
};

TEST(CubatureContructors, ConstructWithoutArray)
{
  // load reference node mesh
  std::unique_ptr<meshBase> mesh1 = meshBase::CreateUnique(nodeMesh);
  // generate gauss point mesh
  std::unique_ptr<GaussCubature> cubeObj = GaussCubature::CreateUnique(mesh1.get()); 
  cubeObj->writeGaussMesh("gaussTestNoData.vtp");
  // load generated mesh
  std::unique_ptr<meshBase> gaussMesh = meshBase::CreateUnique("gaussTestNoData.vtp");
  // load reference gauss mesh
  std::unique_ptr<meshBase> refGaussMesh 
    = meshBase::CreateUnique(refGauss1);
  EXPECT_EQ(0,diffMesh(gaussMesh.get(),refGaussMesh.get()));
} 

// tests construction with array, writeout, loading and interpolation
TEST_F(CubatureTest, InterpolateToGauss)
{
  cuby->writeGaussMesh("gaussTest.vtp");
  std::unique_ptr<meshBase> gaussMesh = meshBase::CreateUnique("gaussTest.vtp");
  std::unique_ptr<meshBase> refGaussMesh = meshBase::CreateUnique(refGauss);
  EXPECT_EQ(0,diffMesh(gaussMesh.get(),refGaussMesh.get()));
}

// test integration methods
TEST_F(CubatureTest, integrateOverAllCells)
{
  std::vector<std::vector<double>> totalIntegralData(cuby->integrateOverAllCells());
  //std::vector<std::vector<double>> refTotalIntegralData(totalIntegralData.size());

  //for (int i = 0; i < totalIntegralData.size(); ++i)
  //{
  //  refTotalIntegralData[i].resize(totalIntegralData[i].size());
  //}
  //refTotalIntegralData[0][0] = 28385.2;
  //refTotalIntegralData[1][0] = 12718.5;
  //refTotalIntegralData[2][0] = 0.00042734;
  //refTotalIntegralData[2][1] = -5.75403e-08;
  //refTotalIntegralData[2][2] = -5.99377e-08;

  std::vector<std::vector<double>> refTotalIntegralData 
    = {{28385.234}, {12718.489}, {10811986}, {0.00042734, -5.754033e-08, -5.99377e-08}};

  for (int i = 0; i < totalIntegralData.size(); ++i)
  {
    for (int j = 0; j < totalIntegralData[i].size(); ++j)
    {
      EXPECT_FLOAT_EQ(totalIntegralData[i][j],refTotalIntegralData[i][j]);
      std::cout << totalIntegralData[i][j] << " ";
    }
    std::cout << std::endl;
  }
  
  std::unique_ptr<meshBase> integralMesh = meshBase::CreateUnique(refInt);
  EXPECT_EQ(0,diffMesh(cuby->getNodeMesh(),integralMesh.get())); 
  
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  nodeMesh = argv[1];
  refGauss = argv[2];
  refGauss1 = argv[3];
	refInt	= argv[4];
  return RUN_ALL_TESTS();
}

//double integrand(const std::vector<double>& coord)
//{
//	return sin(coord[0]) + coord[1]*coord[1]+ cos(coord[2]);
//}
//
//double integrand1(const std::vector<double>& coord)
//{
//	return coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2];
//}
//TEST(IntegrationTest, writeNewCube)
//{
//	std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(cubeMesh);
//	vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
//	da->SetName("integrand");
//	da->SetNumberOfComponents(1);
//	da->SetNumberOfTuples(mesh->getNumberOfPoints());
//
//	for (int i = 0; i < mesh->getNumberOfPoints(); ++i)
//	{
//		double val[1];
//		val[0] = integrand(mesh->getPoint(i));
//		da->InsertTuple(i,val);
//	}
//	mesh->getDataSet()->GetPointData()->AddArray(da);
//	mesh->write();
//}
