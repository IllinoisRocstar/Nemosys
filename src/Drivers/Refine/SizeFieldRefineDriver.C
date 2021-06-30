#include "Drivers/Refine/SizeFieldRefineDriver.H"

#include "meshBase.H"

namespace NEM {
namespace DRV {

SizeFieldRefineDriver::Opts::Opts(Method method, std::string arrayName,
                                  double stdDevMult, bool maxIsMin,
                                  bool transferData)
    : method(method),
      arrayName(std::move(arrayName)),
      stdDevMult(stdDevMult),
      maxIsMin(maxIsMin),
      transferData(transferData) {}

std::string SizeFieldRefineDriver::Opts::getMethodStr() const {
  switch (this->method) {
    case Method::VALUE: return valStr;
    case Method::GRADIENT: return gradStr;
  }
}

SizeFieldRefineDriver::SizeFieldRefineDriver(Files files, Opts opts)
    : RefineDriver(std::move(files)), opts_(std::move(opts)) {}

SizeFieldRefineDriver::SizeFieldRefineDriver()
    : SizeFieldRefineDriver({{}, {}}, {{}, {}, {}, {}, {}}) {}

const SizeFieldRefineDriver::Opts &SizeFieldRefineDriver::getOpts() const {
  return opts_;
}

void SizeFieldRefineDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void SizeFieldRefineDriver::execute() const {
  std::cout << "Size Factor = " << this->opts_.sizeFactor << std::endl;
  std::shared_ptr<meshBase> mesh =
      meshBase::CreateShared(this->files_.inputMeshFile);
  std::cout << "\n";
  mesh->report();
  std::cout << "\n";
  // Edge scale is unused
  mesh->refineMesh(this->opts_.getMethodStr(), this->opts_.arrayName,
                   this->opts_.stdDevMult, this->opts_.maxIsMin, {},
                   this->files_.outputMeshFile, this->opts_.transferData,
                   this->opts_.sizeFactor);
}

}  // namespace DRV
}  // namespace NEM
