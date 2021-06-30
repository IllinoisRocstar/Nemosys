#include "Drivers/Refine/UniformRefineDriver.H"

#include "meshBase.H"

namespace NEM {
namespace DRV {

UniformRefineDriver::Opts::Opts(bool transferData, double edgeScale)
    : transferData(transferData), edgeScale(edgeScale) {}

UniformRefineDriver::UniformRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(opts) {}

UniformRefineDriver::UniformRefineDriver()
    : UniformRefineDriver({{}, {}}, {{}, {}}) {}

const UniformRefineDriver::Opts &UniformRefineDriver::getOpts() const {
  return opts_;
}

void UniformRefineDriver::setOpts(Opts opts) { this->opts_ = opts; }

void UniformRefineDriver::execute() const {
  std::shared_ptr<meshBase> mesh =
      meshBase::CreateShared(this->files_.inputMeshFile);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  mesh->refineMesh(Opts::method, this->opts_.edgeScale,
                   this->files_.outputMeshFile, this->opts_.transferData);
}

}  // namespace DRV
}  // namespace NEM
