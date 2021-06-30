#include "Drivers/MeshGen/CFMeshMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "cfmeshGen.H"
#include "meshBase.H"

namespace NEM {
namespace DRV {

CFMeshMeshGenDriver::Opts::Opts(cfmeshParams params)
    : params(std::move(params)) {}

CFMeshMeshGenDriver::CFMeshMeshGenDriver(Files files, cfmeshParams params)
    : files_(std::move(files)), opts_(std::move(params)) {}

CFMeshMeshGenDriver::CFMeshMeshGenDriver()
    : CFMeshMeshGenDriver({{}, {}}, {}) {}

const CFMeshMeshGenDriver::Files &CFMeshMeshGenDriver::getFiles() const {
  return files_;
}

void CFMeshMeshGenDriver::setFiles(Files files) {
  this->opts_.params.geomFilePath = files.inputGeoFile;
  this->files_ = std::move(files);
}

const cfmeshParams &CFMeshMeshGenDriver::getParams() const {
  return getOpts().params;
}

void CFMeshMeshGenDriver::setParams(cfmeshParams params) {
  setOpts(Opts{std::move(params)});
}

const CFMeshMeshGenDriver::Opts &CFMeshMeshGenDriver::getOpts() const {
  return opts_;
}

void CFMeshMeshGenDriver::setOpts(Opts opts) {
  if (!opts.params.geomFilePath.empty()) {
    this->files_.inputGeoFile = opts.params.geomFilePath;
    this->opts_ = std::move(opts);
  } else {
    this->opts_ = std::move(opts);
    this->opts_.params.geomFilePath = this->files_.inputGeoFile;
  }
}

void CFMeshMeshGenDriver::execute() const {
  auto paramsCopy = this->opts_.params;
  cfmeshGen generator{&paramsCopy};
  // TODO: Make sure blockMeshGen::createMeshFromSTL sets return value and check
  //  it here
  // Parameter not used
  generator.createMeshFromSTL(nullptr);
  std::string newname = nemAux::trim_fname(this->files_.inputGeoFile, ".vtu");
  auto mesh =
      meshBase::CreateShared(meshBase::Create(generator.getDataSet(), newname));
  mesh->setFileName(this->files_.outputMeshFile);
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
