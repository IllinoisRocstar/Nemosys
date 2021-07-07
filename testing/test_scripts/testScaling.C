#include <gtest/gtest.h>
#include <Integration/Cubature.H>
#include <Drivers/TransferDriver.H>
#include <SolutionVerification/OrderOfAccuracy.H>
#include <Mesh/meshBase.H>

#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include <omp.h>

const char *coarseFname;
const char *fineFname;
const char *finerFname;

std::shared_ptr<meshBase> coarse;
std::shared_ptr<meshBase> fine;
std::shared_ptr<meshBase> finer;

GaussCubature *cubature;
OrderOfAccuracy *oac;

TEST(ScalingTest, Read) {
  omp_set_num_threads(4);
  std::cout << "MAX NUM THREADS : " << omp_get_max_threads() << std::endl;
  std::cout << "THREAD NUM : " << omp_get_thread_num() << std::endl;
  coarse = meshBase::CreateShared(coarseFname);
  fine = meshBase::CreateShared(fineFname);
  finer = meshBase::CreateShared(finerFname);
}

TEST(ScalingTest, SetCoarseData) {
  // add linear field to coarse grid
  auto linearField = vtkSmartPointer<vtkDoubleArray>::New();
  linearField->SetName("linear field");
  for (int i = 0; i < coarse->getDataSet()->GetNumberOfPoints(); ++i) {
    double x[3];
    coarse->getDataSet()->GetPoint(i, x);
    linearField->InsertValue(i, x[0]);
  }
  coarse->getDataSet()->GetPointData()->AddArray(linearField);

  /*
  // add quadratic field to coarse grid
  auto quadraticField = vtkSmartPointer<vtkDoubleArray>::New();
  quadraticField->SetName("quadratic field");
  for(int i = 0; i < coarse->getDataSet()->GetNumberOfPoints(); ++i) {
    double x[3];
    coarse->getDataSet()->GetPoint(i, x);
    quadraticField->InsertValue(i, x[0]*x[0]);
  }
  coarse->getDataSet()->GetPointData()->AddArray(quadraticField);
  */
}

/*
TEST(ScalingTest, SetFineData) {
  // add linear field to fine grid
  auto linearField = vtkSmartPointer<vtkDoubleArray>::New();
  linearField->SetName("linear field");
  for (int i = 0; i < fine->getDataSet()->GetNumberOfPoints(); ++i) {
    double x[3];
    fine->getDataSet()->GetPoint(i, x);
    linearField->InsertValue(i, x[0]);
  }
  fine->getDataSet()->GetPointData()->AddArray(linearField);

  // add quadratic field to fine grid
  auto quadraticField = vtkSmartPointer<vtkDoubleArray>::New();
  quadraticField->SetName("quadratic field");
  for (int i = 0; i < fine->getDataSet()->GetNumberOfPoints(); ++i) {
    double x[3];
    fine->getDataSet()->GetPoint(i, x);
    quadraticField->InsertValue(i, x[0] * x[0]);
  }
  fine->getDataSet()->GetPointData()->AddArray(quadraticField);
}

TEST(ScalingTest, SetFinerData) {
  // add linear field to finer grid
  auto linearField = vtkSmartPointer<vtkDoubleArray>::New();
  linearField->SetName("linear field");
  for (int i = 0; i < finer->getDataSet()->GetNumberOfPoints(); ++i) {
    double x[3];
    finer->getDataSet()->GetPoint(i, x);
    linearField->InsertValue(i, x[0]);
  }
  finer->getDataSet()->GetPointData()->AddArray(linearField);

  // add quadratic field to finer grid
  auto quadraticField = vtkSmartPointer<vtkDoubleArray>::New();
  quadraticField->SetName("quadratic field");
  for (int i = 0; i < finer->getDataSet()->GetNumberOfPoints(); ++i) {
    double x[3];
    finer->getDataSet()->GetPoint(i, x);
    quadraticField->InsertValue(i, x[0] * x[0]);
  }
  finer->getDataSet()->GetPointData()->AddArray(quadraticField);
}
*/

TEST(ScalingTest, Transfer) {
  // transfer
  std::string method("Consistent Interpolation");
  auto transfer = NEM::DRV::TransferDriver::CreateTransferObject(
      coarse.get(), fine.get(), method);
  std::cerr << "transferring..." << std::endl;
  transfer->run();

  // check target values
  auto transferredValues = vtkDoubleArray::SafeDownCast(
      fine->getDataSet()->GetPointData()->GetArray(0));
  for (int i = 0; i < fine->getDataSet()->GetNumberOfPoints(); ++i) {
    double val = transferredValues->GetValue(i);
    double x[3];
    fine->getDataSet()->GetPoint(i, x);
    ASSERT_NEAR(x[0], val, 1e-14);
  }
}

/*
TEST(ScalingTest, CubatureConstruction) {
  const std::vector<int> arrayIDs = {0};
  cubature = new GaussCubature(finer.get(), arrayIDs);
}

TEST(ScalingTest, IntegrateOverMesh) {
  const std::vector<int> arrayIDs = {0};
  finer->integrateOverMesh(arrayIDs);
  // std::vector<std::vector<double>>
  integral(coarse->integrateOverMesh(arrayIDs));
}

TEST(ScalingTest, VerificationConstruction) {
  const std::vector<int> arrayIDs = {0};
  oac = new OrderOfAccuracy(finer.get(), fine.get(), coarse.get(), arrayIDs);
}

TEST(ScalingTest, CheckAysmptoticRange) {
  std::vector<std::vector<double>> ratios = oac->checkAsymptoticRange();
}
*/

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 4);
  coarseFname = argv[1];
  fineFname = argv[2];
  finerFname = argv[3];
  return RUN_ALL_TESTS();
}
