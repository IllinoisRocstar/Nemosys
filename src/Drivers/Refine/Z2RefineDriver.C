#include "Drivers/Refine/Z2RefineDriver.H"

#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

Z2RefineDriver::Opts::Opts(bool transferData, std::string arrayName, int order)
    : transferData(transferData),
      arrayName(std::move(arrayName)),
      order(order) {}

Z2RefineDriver::Z2RefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

Z2RefineDriver::Z2RefineDriver() : Z2RefineDriver({{}, {}}, {{}, {}, {}}) {}

const Z2RefineDriver::Opts &Z2RefineDriver::getOpts() const {
  return opts_;
}

void Z2RefineDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void Z2RefineDriver::execute() const {
  std::shared_ptr<meshBase> mesh =
      meshBase::CreateShared(this->files_.inputMeshFile);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(Opts::method, this->opts_.arrayName, this->opts_.order,
                   this->files_.outputMeshFile, this->opts_.transferData);
}

}  // namespace DRV
}  // namespace NEM
