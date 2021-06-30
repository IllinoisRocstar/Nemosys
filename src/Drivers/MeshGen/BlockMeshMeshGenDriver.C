#include "Drivers/MeshGen/BlockMeshMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "blockMeshGen.H"
#include "meshBase.H"

namespace NEM {
namespace DRV {

BlockMeshMeshGenDriver::Opts::Opts(blockMeshParams params)
    : params(std::move(params)) {}

BlockMeshMeshGenDriver::BlockMeshMeshGenDriver(Files file,
                                               blockMeshParams params)
    : file_(std::move(file)), opts_(std::move(params)) {}

BlockMeshMeshGenDriver::BlockMeshMeshGenDriver()
    : BlockMeshMeshGenDriver(Files{{}}, {}) {}

const BlockMeshMeshGenDriver::Files &BlockMeshMeshGenDriver::getFiles() const {
  return file_;
}

void BlockMeshMeshGenDriver::setFiles(Files file) {
  this->file_ = std::move(file);
}

const blockMeshParams &BlockMeshMeshGenDriver::getParams() const {
  return getOpts().params;
}

void BlockMeshMeshGenDriver::setParams(blockMeshParams params) {
  setOpts(Opts{std::move(params)});
}

const BlockMeshMeshGenDriver::Opts &BlockMeshMeshGenDriver::getOpts() const {
  return opts_;
}

void BlockMeshMeshGenDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void BlockMeshMeshGenDriver::execute() const {
  auto paramsCopy = this->opts_.params;
  blockMeshGen generator{&paramsCopy};
  // TODO: Make sure blockMeshGen::createMeshFromSTL sets return value and check
  //  it here
  // Parameter not used
  generator.createMeshFromSTL(nullptr);
  auto mesh = meshBase::CreateShared(
      meshBase::Create(generator.getDataSet(), this->file_.outputFile));
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
