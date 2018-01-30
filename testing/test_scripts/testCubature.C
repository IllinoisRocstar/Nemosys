#include <Cubature.H>
#include <gtest.h>

int diffVTU(meshBase* mesh1, meshBase* mesh2)
{
  double tol = 1e-6;

  if (mesh1->getNumberOfPoints() != mesh2->getNumberOfPoints() ||
      mesh1->getNumberOfCells() != mesh2->getNumberOfCells())
  {
    std::cerr << "Meshes don't have the same number of points or cells" << std::endl;
    exit(1);
  }

  for (int i = 0; i < mesh1->getNumberOfPoints(); ++i)
  {
    std::vector<double> coord1 = mesh1->getPoint(i);
    std::vector<double> coord2 = mesh2->getPoint(i);
    for (int j = 0; j < 3; ++j)
    {
      if (std::fabs(coord1[j]-coord2[j]) > tol)
      {
        std::cerr << "Meshes differ in point coordinates" << std::endl;
        exit(1);
      }
    } 
  }

  for (int i = 0; i < mesh1->getNumberOfCells(); ++i)
  {
    std::vector<std::vector<double>> cell1 = mesh1->getCellVec(i);
    std::vector<std::vector<double>> cell2 = mesh2->getCellVec(i);
    if (cell1.size() != cell2.size())
    {
      std::cerr << "Meshes differ in cells" << std::endl;
      exit(1); 
    }
    for (int j = 0; j < cell1.size(); ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        if (std::fabs(cell1[j][k] - cell2[j][k]) > tol)
        {
          std::cerr << "Meshes differ in cells" << std::endl;
          exit(1); 
        }
      }
    }   
  }

  vtkSmartPointer<vtkPointData> pd1 = vtkSmartPointer<vtkPointData>::New();
  pd1 = mesh1->getDataSet()->GetPointData();
  vtkSmartPointer<vtkPointData> pd2 = vtkSmartPointer<vtkPointData>::New();
  pd2 = mesh2->getDataSet()->GetPointData(); 
  int numArr1 = pd1->GetNumberOfArrays(); 
  int numArr2 = pd2->GetNumberOfArrays(); 
 
  if (numArr1 != numArr2)
  {
    std::cerr << "Meshes have different numbers of point data" << std::endl;
    exit(1);
  }

  for (int i = 0; i < numArr1; ++i)
  {
    vtkDataArray* da1 = pd1->GetArray(i);
    vtkDataArray* da2 = pd2->GetArray(i);
    int numComponent = da1->GetNumberOfComponents();
    for (int j = 0; j < mesh1->getNumberOfPoints(); ++j)
    {
      double comps1[numComponent];
      double comps2[numComponent];
      da1->GetTuple(j, comps1);
      da2->GetTuple(j, comps2);
      for (int k = 0; k < numComponent; ++k)
      {
        if (std::fabs(comps1[k] - comps2[k]) > tol)
        {
          std::cerr << "Meshes differ in point data values" << std::endl;
          exit(1);
        }
      }
    }
  }
  if (mesh1) delete mesh1;
  if (mesh2) delete mesh2;
  std::cerr << "Meshes are the same" << std::endl;
  return EXIT_SUCCESS;
}

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
      GaussCubature* cuby = new  GaussCubature(mesh);
    }
  
    virtual ~CubatureTest()
    {
      delete mesh;
      delete cuby;
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

TEST_F(CubatureTest, InterpolateToGauss)
{
  cuby->interpolateToGaussPoints(arrayIDs);   
  cuby->writeGaussMesh();
  meshBase* gaussMesh = meshBase::Create("gaussTest.vtp");
  meshBase* refGaussMesh = meshBase::Create(refGauss);
  EXPECT_EXIT(diffVTU(gaussMesh,refGaussMesh), ::testing::ExitedWithCode(0),
              "Meshes are the same"); 
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 3);
  nodeMesh = argv[1];
  refGauss = argv[2];
  return RUN_ALL_TESTS();
}

