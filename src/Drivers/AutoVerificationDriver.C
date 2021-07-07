#include "Drivers/AutoVerificationDriver.H"

#include "SolutionVerification/OrderOfAccuracy.H"
#include "Mesh/meshBase.H"

#ifdef HAVE_OPENMP
#  include <omp.h>
#endif

namespace NEM {
namespace DRV {

AutoVerificationDriver::Files::Files(std::string coarseMesh,
                                     std::string fineMesh,
                                     std::string finerMesh)
    : coarseMeshFile(std::move(coarseMesh)),
      fineMeshFile(std::move(fineMesh)),
      finerMeshFile(std::move(finerMesh)) {}

AutoVerificationDriver::Opts::Opts(std::vector<int> arrayIds)
    : arrayIds(std::move(arrayIds)) {}

AutoVerificationDriver::AutoVerificationDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

AutoVerificationDriver::AutoVerificationDriver()
    : AutoVerificationDriver({{}, {}, {}}, Opts{{}}) {}

const AutoVerificationDriver::Files &AutoVerificationDriver::getFiles() const {
  return files_;
}

void AutoVerificationDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const AutoVerificationDriver::Opts &AutoVerificationDriver::getOpts() const {
  return opts_;
}

void AutoVerificationDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

jsoncons::string_view AutoVerificationDriver::getProgramType() const {
  return programType;
}

void AutoVerificationDriver::execute() const {
  auto coarseMesh = meshBase::CreateShared(this->files_.coarseMeshFile);
  auto fineMesh = meshBase::CreateShared(this->files_.fineMeshFile);
  auto finerMesh = meshBase::CreateShared(this->files_.finerMeshFile);
#ifdef HAVE_OPENMP
  auto threads = this->opts_.numThreads.value_or(omp_get_max_threads());
  omp_set_num_threads(threads);
  std::cout << "Number of threads set to : " << threads << std::endl;
#else
  if (this->opts_.numThreads) {
    std::cerr << "OpenMP is not enabled. Verification will continue in serial."
              << std::endl;
  }
#endif
  std::cout << "Running verification." << std::endl;
  auto oac = std::make_shared<OrderOfAccuracy>(
      coarseMesh.get(), fineMesh.get(), finerMesh.get(), this->opts_.arrayIds,
      this->opts_.transferType, this->opts_.targetGCI);
  std::cout << "Checking if in asymptotic range." << std::endl;
  std::cout << "Target GCI is set to : " << oac->getTargetGCI() << std::endl;
  auto asymp = oac->checkAsymptoticRange();
  bool inRange = true;
  for (int i = 0; i < asymp.size(); ++i) {
    for (int j = 0; j < asymp[0].size(); ++j) {
      double gci = asymp[i][j];
      if (gci > oac->getTargetGCI()) {
        std::cout << "GCI of " << gci;
        std::cout << " exceeds target GCI of " << oac->getTargetGCI()
                  << std::endl;
        std::cout << "at array " << i << ", component " << j << std::endl;
        inRange = false;
      }
    }
  }
  if (inRange) {
    std::cout << "Grid is in target asymptotic range." << std::endl;
  } else {
    std::cout << "Grid is not in target asymptotic range." << std::endl;
  }
}

}  // namespace DRV
}  // namespace NEM
