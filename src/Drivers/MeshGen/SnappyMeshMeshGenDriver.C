#include "Drivers/MeshGen/SnappyMeshMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "Mesh/meshBase.H"
#include "MeshGeneration/snappymeshGen.H"

namespace NEM {
namespace DRV {

SnappyMeshMeshGenDriver::Opts::Opts(snappymeshParams params)
    : params(std::move(params)) {}

SnappyMeshMeshGenDriver::SnappyMeshMeshGenDriver(Files files,
                                                 snappymeshParams params)
    : files_(std::move(files)), opts_(std::move(params)) {}

SnappyMeshMeshGenDriver::SnappyMeshMeshGenDriver()
    : SnappyMeshMeshGenDriver({{}, {}}, {}) {}

const SnappyMeshMeshGenDriver::Files &SnappyMeshMeshGenDriver::getFiles()
    const {
  return files_;
}

void SnappyMeshMeshGenDriver::setFiles(Files files) {
  this->opts_.params.geomFileName = files.inputGeoFile;
  this->files_ = std::move(files);
}

const snappymeshParams &SnappyMeshMeshGenDriver::getParams() const {
  return getOpts().params;
}

void SnappyMeshMeshGenDriver::setParams(snappymeshParams params) {
  setOpts(Opts{std::move(params)});
}

const SnappyMeshMeshGenDriver::Opts &SnappyMeshMeshGenDriver::getOpts() const {
  return opts_;
}

void SnappyMeshMeshGenDriver::setOpts(Opts opts) {
  if (!opts.params.geomFileName.empty()) {
    this->files_.inputGeoFile = opts.params.geomFileName;
    this->opts_ = std::move(opts);
  } else {
    this->opts_ = std::move(opts);
    this->opts_.params.geomFileName = this->files_.inputGeoFile;
  }
}

void SnappyMeshMeshGenDriver::execute() const {
  auto paramsCopy = this->opts_.params;
  snappymeshGen generator{&paramsCopy};
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
