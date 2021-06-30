#include "Drivers/Refine/FoamRefineDriver.H"

#include <pointFieldsFwd.H>

#include <interpolatePointToCell.H>
#include "AMRFoam.H"
#include "AuxiliaryFunctions.H"
#include "MeshManipulationFoam.H"
#include "foamMesh.H"
#include "meshBase.H"
#include "vtkMesh.H"

#ifdef MLAMR
#  include <fdeep/fdeep.hpp>
#endif

// helpers for execute
namespace {
std::unique_ptr<Foam::AMRFoam> initializeRefine(
    const NEM::DRV::RefineDriver::Files &files,
    const NEM::DRV::FoamRefineOptsBase &opts, const Foam::Time &runTime,
    std::shared_ptr<meshBase> mesh = nullptr) {
  if (mesh == nullptr) {
    mesh = meshBase::CreateShared(files.inputMeshFile);
  }

  // Reading mesh from vtk file.
  FOAM::foamMesh *fm = new FOAM::foamMesh(mesh);
  fm->write(files.outputMeshFile);

  Foam::polyMesh mesh1(Foam::IOobject(Foam::polyMesh::defaultRegion,
                                      runTime.timeName(), runTime,
                                      Foam::IOobject::MUST_READ));

  // Initialized AMR class
  auto amr = std::unique_ptr<Foam::AMRFoam>{new Foam::AMRFoam(mesh1)};

  if (opts.writeFieldData)
    if (!opts.writeMesh) amr->enableUpdatedField();

  if (opts.writeRefHistory)
    if (!opts.writeMesh) amr->enableRefHistoryData();

  if (opts.writeMesh) amr->enableMeshWriting();

  return amr;
}

void finishRefine(Foam::AMRFoam *amr, const std::string &outputMeshFile,
                  const Foam::Time &runTime) {
  amr->writeMesh();

  // Export final mesh to vtu format
  FOAM::foamMesh *fm2 = new FOAM::foamMesh(false);
  fm2->readAMR(runTime);
  vtkMesh *vm = new vtkMesh(fm2->getDataSet(), outputMeshFile);
  vm->report();
  vm->write();
  delete vm;
}
}  // namespace

namespace NEM {
namespace DRV {

FoamRefineDriver::Opts::Opts(Criteria refCriteria, std::string inputFieldFile)
    : FoamRefineOptsBase(),
      inputFieldFile(std::move(inputFieldFile)),
      refCriteria(refCriteria) {}

FoamRefineDriver::FoamRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

FoamRefineDriver::FoamRefineDriver() : FoamRefineDriver({{}, {}}, {{}, {}}) {}

const FoamRefineDriver::Opts &FoamRefineDriver::getOpts() const {
  return opts_;
}

void FoamRefineDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

void FoamRefineDriver::execute() const {
  // Initializes AMR workflow
  auto *mshManip = new MeshManipulationFoam();
  mshManip->initAMRWorkFlow(this->opts_.startTime, this->opts_.timeStep,
                            this->opts_.endTime);
  if (mshManip) delete mshManip;

  // Starting
  int argc = 1;
  char **argv = new char *[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info << "Create time\n" << Foam::endl;
  Foam::Time runTime(Foam::Time::controlDictName, args);

  auto amr = initializeRefine(this->files_, this->opts_, runTime);

  int nSteps = opts_.endTime / opts_.timeStep;

  // Creates scalar field from incoming text/CSV files
  volScalarField meshFieldXY =
      amr->readIncomingCellField(this->opts_.inputFieldFile, "meshFieldXY");

  if (opts_.refCriteria == Opts::Criteria::VALUE) {
    for (int i = 0; i < nSteps; i++) {
      amr->updateAMR(this->opts_.refineInterval, this->opts_.maxRefinement,
                     meshFieldXY, this->opts_.lowerRefineLevel,
                     this->opts_.upperRefineLevel, this->opts_.unrefineAbove,
                     this->opts_.unrefineBelow, this->opts_.nBufferLayers,
                     this->opts_.maxCells);

      runTime++;
    }
    meshFieldXY.write();
  } else {  // opts_.refCriteria == Opts::Criteria::GRADIENT
    volScalarField magGrad = amr->getGradient(meshFieldXY);
    for (int i = 0; i < nSteps; i++) {
      amr->updateAMR(this->opts_.refineInterval, this->opts_.maxRefinement,
                     magGrad, this->opts_.lowerRefineLevel,
                     this->opts_.upperRefineLevel, this->opts_.unrefineAbove,
                     this->opts_.unrefineBelow, this->opts_.nBufferLayers,
                     this->opts_.maxCells);

      runTime++;
    }
    magGrad.write();
  }

  finishRefine(amr.get(), this->files_.outputMeshFile, runTime);

  std::cout << "End!" << std::endl;
}

#ifdef MLAMR
FoamMLRefineDriver::FoamMLRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

FoamMLRefineDriver::FoamMLRefineDriver() : FoamMLRefineDriver({{}, {}}, {}) {}

const FoamMLRefineDriver::Opts &FoamMLRefineDriver::getOpts() const {
  return opts_;
}

void FoamMLRefineDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void FoamMLRefineDriver::execute() const {
  // Initializes AMR workflow
  auto *mshManip = new MeshManipulationFoam();
  mshManip->initAMRWorkFlow(this->opts_.startTime, this->opts_.timeStep,
                            this->opts_.endTime);
  if (mshManip) delete mshManip;

  // Starting
  int argc = 1;
  char **argv = new char *[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info << "Create time\n" << Foam::endl;
  Foam::Time runTime(Foam::Time::controlDictName, args);

  std::shared_ptr<meshBase> mesh =
      meshBase::CreateShared(this->files_.inputMeshFile);
  auto amr = initializeRefine(this->files_, this->opts_, runTime, mesh);

  int nSteps = opts_.endTime / opts_.timeStep;

  std::vector<double> nonDimUGrad;
  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> Z;
  mesh->getCellDataArray("nonDimUGrad", nonDimUGrad);
  mesh->getCellDataArray("X", X);
  mesh->getCellDataArray("Y", Y);
  mesh->getCellDataArray("Z", Z);

  // Loading ML Model
  const auto model = fdeep::load_model(this->opts_.mlModel);

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
  volScalarField meshField = amr->assignToVolScalarField(refinementVec);

  for (int i = 0; i < nSteps; i++) {
    amr->updateAMRML(this->opts_.refineInterval, this->opts_.maxRefinement,
                     this->opts_.nBufferLayers, this->opts_.maxCells,
                     meshField);

    runTime++;
  }

  finishRefine(amr.get(), this->files_.outputMeshFile, runTime);
  std::cout << "End!" << std::endl;
}
#endif

}  // namespace DRV
}  // namespace NEM
