#include <Cubature.H>
#include <gtest.h>

const char* nodeMesh;
const char* refGauss;
const char* refGauss1;
const char* cubeMesh;

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

  std::vector<std::vector<double>> totalIntegralData(cuby->integrateOverAllCells());
  std::string fname("integrationTest.vtu");
  cuby->getNodeMesh()->write(fname);

  for (int i = 0; i < totalIntegralData.size(); ++i)
  {
    for (int j = 0; j < totalIntegralData[i].size(); ++j)
    {
      std::cout << totalIntegralData[i][j] << " ";
    }
    std::cout << std::endl;
  }

}

double integrand(const std::vector<double>& coord)
{
	return sin(coord[0]) + coord[1]*coord[1]+ cos(coord[2]);
}

double integrand1(const std::vector<double>& coord)
{
	return coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2];
}

TEST(IntegrationTest, writeNewCube)
{
	std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(cubeMesh);
	vtkSmartPointer<vtkDoubleArray> da = vtkSmartPointer<vtkDoubleArray>::New();
	da->SetName("integrand");
	da->SetNumberOfComponents(1);
	da->SetNumberOfTuples(mesh->getNumberOfPoints());

	for (int i = 0; i < mesh->getNumberOfPoints(); ++i)
	{
		double val[1];
		val[0] = integrand(mesh->getPoint(i));
		da->InsertTuple(i,val);
	}
	mesh->getDataSet()->GetPointData()->AddArray(da);
	mesh->write("sphereMeshWithIntegrand.vtu");
}
// test integration methods
TEST(IntegrationTest, integrateOverCube)
{
	std::unique_ptr<meshBase> mesh = meshBase::CreateUnique(cubeMesh);
	std::vector<int> arrayIDs = {0,1};
	std::unique_ptr<GaussCubature> cuby = GaussCubature::CreateUnique(mesh.get(),arrayIDs);
	std::vector<std::vector<double>> integral(cuby->integrateOverAllCells());
	std::cout << "Integral value 1: " << integral[0][0] << std::endl;
	std::cout << "Integral value 2: " << integral[1][0] << std::endl;
	mesh->write("sphereMeshWithIntegral.vtu");
}



int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  nodeMesh = argv[1];
  refGauss = argv[2];
  refGauss1 = argv[3];
	cubeMesh	= argv[4];
  return RUN_ALL_TESTS();
}

