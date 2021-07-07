#include "Drivers/MeshGen/NetgenMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "Mesh/meshBase.H"
#include "MeshGeneration/netgenGen.H"

namespace NEM {
namespace DRV {

NetgenMeshGenDriver::Opts::Opts(netgenParams params)
    : params(std::move(params)) {}

NetgenMeshGenDriver::NetgenMeshGenDriver(Files files, netgenParams params)
    : files_(std::move(files)), opts_(std::move(params)) {}

NetgenMeshGenDriver::NetgenMeshGenDriver()
    : NetgenMeshGenDriver({{}, {}}, {}) {}

const NetgenMeshGenDriver::Files &NetgenMeshGenDriver::getFiles() const {
  return files_;
}

void NetgenMeshGenDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const netgenParams &NetgenMeshGenDriver::getParams() const {
  return getOpts().params;
}

void NetgenMeshGenDriver::setParams(netgenParams params) {
  setOpts(Opts{std::move(params)});
}

const NetgenMeshGenDriver::Opts &NetgenMeshGenDriver::getOpts() const {
  return opts_;
}

void NetgenMeshGenDriver::setOpts(Opts opts) { this->opts_ = std::move(opts); }

void NetgenMeshGenDriver::execute() const {
  netgenGen generator{&this->opts_.params};
  int status = generator.createMeshFromSTL(this->files_.inputGeoFile.c_str());
  if (status) {
    std::cerr << "Mesh Generation encountered error." << std::endl;
    exit(1);
  }
  std::string newname = nemAux::trim_fname(this->files_.inputGeoFile, ".vol");
  auto mesh = meshBase::CreateShared(meshBase::exportVolToVtk(newname));
  mesh->setFileName(this->files_.outputMeshFile);
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
