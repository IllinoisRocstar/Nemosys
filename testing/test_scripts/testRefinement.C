#include <RefineDriver.H>
#include <gtest.h>

#include "meshBase.H"

#ifdef HAVE_CFMSH
#  include "AMRFoam.H"
#endif
#ifdef MLAMR
#  include <fdeep/fdeep.hpp>
#  include "vtkMesh.H"
#endif

const char *refineValueJSON;
const char *refineValueVTU;
const char *refineValueGoldVTU;
const char *refineUniformJSON;
const char *refineUniformVTU;
const char *refineUniformGoldVTU;
const char *AMRValueJSON;
const char *mlPredictVTU;
const char *mlModelFile;

int genTest(const char *jsonF, const char *newName, const char *goldName) {
  std::string fname(jsonF);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  std::unique_ptr<NEM::DRV::RefineDriver> refineDriver =
      std::unique_ptr<NEM::DRV::RefineDriver>(
          NEM::DRV::RefineDriver::readJSON(inputjson));

  std::unique_ptr<meshBase> goldMesh = meshBase::CreateUnique(goldName);
  std::unique_ptr<meshBase> newMesh = meshBase::CreateUnique(newName);

  // return 0;
  return diffMesh(goldMesh.get(), newMesh.get());
}

#ifdef HAVE_CFMSH
int AMRTest(const char *jsonF2) {
  std::string fname(jsonF2);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << jsonF2 << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;
  std::unique_ptr<NEM::DRV::RefineDriver> refineDriver =
      std::unique_ptr<NEM::DRV::RefineDriver>(
          NEM::DRV::RefineDriver::readJSON(inputjson));

  std::unique_ptr<meshBase> genMesh = meshBase::CreateUnique("refined_AMR.vtu");
  std::unique_ptr<meshBase> refMesh =
      meshBase::CreateUnique("refined_AMR_ref.vtu");

  if ((genMesh->getNumberOfPoints() >= refMesh->getNumberOfPoints() * 0.96) &&
      genMesh->getNumberOfPoints() <= refMesh->getNumberOfPoints() * 1.04)
    if ((genMesh->getNumberOfCells() >= refMesh->getNumberOfCells() * 0.96) &&
        genMesh->getNumberOfCells() <= refMesh->getNumberOfCells() * 1.04)
      return 0;
    else
      return 1;
  else
    return 1;
}
#endif

#ifdef MLAMR
// ML Test
int MLAMRTest(const std::string MLName, const std::string MeshName) {
  std::string fname(MLName);
  std::ifstream inputStream(fname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << MLName << std::endl;
    exit(1);
  }

  // Read incoming vtk mesh
  std::shared_ptr<meshBase> mb = meshBase::CreateShared(MeshName);
  vtkMesh *vm = new vtkMesh(mb->getDataSet(), "outMesh.vtu");

  // Get ML inputs and reference output
  std::vector<double> nonDimUGrad;
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  std::vector<double> refOutput;

  vm->getCellDataArray("nonDimUGrad", nonDimUGrad);
  vm->getCellDataArray("X", X);
  vm->getCellDataArray("Y", Y);
  vm->getCellDataArray("Z", Z);
  vm->getCellDataArray("ReferenceOutput", refOutput);

  // Loading ML Model
  const auto model = fdeep::load_model(MLName);

  // Making predictions
  std::vector<int> refinementVec;
  for (int i = 0; i < nonDimUGrad.size(); i++) {
    refinementVec.push_back(0);
  }
  for (int i = 0; i < nonDimUGrad.size(); i++) {
    const auto result = model.predict(
        {fdeep::tensor(fdeep::tensor_shape(static_cast<double>(4)),
                       {X[i], Y[i], Z[i], nonDimUGrad[i]})});

    if (result[0].get(0, 0, 0, 0, 0) >= 0.5) refinementVec[i] = 1;
  }

  bool compare = true;

  for (int i = 0; i < refinementVec.size(); i++) {
    if (refinementVec[i] != (int)refOutput[i]) compare = false;
  }

  if (compare)
    return 0;
  else
    return 1;
}
#endif

TEST(RefinementDriverTest, RefineValue) {
  EXPECT_EQ(0, genTest(refineValueJSON, refineValueVTU, refineValueGoldVTU));
}

TEST(RefinementDriverTest, RefineUniform) {
  EXPECT_EQ(0,
            genTest(refineUniformJSON, refineUniformVTU, refineUniformGoldVTU));
}

#ifdef HAVE_CFMSH
TEST(AdaptiveMeshRefinementTest, Refine_Unrefine) {
  EXPECT_EQ(0, AMRTest(AMRValueJSON));
}
#endif

#ifdef MLAMR
TEST(MachineLearningTest, ML_Predictions) {
  EXPECT_EQ(0, MLAMRTest(mlModelFile, mlPredictVTU));
}
#endif

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 10);
  refineValueJSON = argv[1];
  refineValueVTU = argv[2];
  refineValueGoldVTU = argv[3];
  refineUniformJSON = argv[4];
  refineUniformVTU = argv[5];
  refineUniformGoldVTU = argv[6];
  AMRValueJSON = argv[7];
  mlPredictVTU = argv[8];
  mlModelFile = argv[9];
  return RUN_ALL_TESTS();
}
