#include "Drivers/MeshGen/GmshMeshGenDriver.H"

#include "AuxiliaryFunctions.H"
#include "MeshGeneration/gmshGen.H"
#include "Mesh/meshBase.H"

namespace NEM {
namespace DRV {

GmshMeshGenDriver::Opts::Opts(NEM::GEN::gmshParams params)
    : params(std::move(params)) {}

GmshMeshGenDriver::GmshMeshGenDriver(Files files, NEM::GEN::gmshParams params)
    : files_(std::move(files)), opts_(std::move(params)) {}

GmshMeshGenDriver::GmshMeshGenDriver() : GmshMeshGenDriver({{}, {}}, {}) {}

const GmshMeshGenDriver::Files &GmshMeshGenDriver::getFiles() const {
  return files_;
}

void GmshMeshGenDriver::setFiles(Files files) {
  this->opts_.params.ofname = files.outputMeshFile;
  this->files_ = std::move(files);
}

const NEM::GEN::gmshParams &GmshMeshGenDriver::getParams() const {
  return getOpts().params;
}

void GmshMeshGenDriver::setParams(NEM::GEN::gmshParams params) {
  setOpts(Opts{std::move(params)});
}

const GmshMeshGenDriver::Opts &GmshMeshGenDriver::getOpts() const {
  return opts_;
}

void GmshMeshGenDriver::setOpts(Opts opts) {
  if (!opts.params.ofname.empty()) {
    this->files_.outputMeshFile = opts.params.ofname;
    this->opts_ = std::move(opts);
  } else {
    this->opts_ = std::move(opts);
    this->opts_.params.ofname = this->files_.outputMeshFile;
  }
}

void GmshMeshGenDriver::execute() const {
  auto paramsCopy = this->opts_.params;
  NEM::GEN::gmshGen generator{&paramsCopy};
  int status = generator.createMeshFromSTL(this->files_.inputGeoFile.c_str());
  if (status) {
    std::cerr << "Mesh Generation encountered error." << std::endl;
    exit(1);
  }
  std::string outputType = nemAux::find_ext(this->files_.outputMeshFile);
  if (outputType == ".msh") {
    return;
  }
  std::string newname = nemAux::trim_fname(this->files_.inputGeoFile, ".msh");
  auto mesh = meshBase::CreateShared(meshBase::exportGmshToVtk(newname));
  mesh->setFileName(this->files_.outputMeshFile);
  mesh->report();
  mesh->write();
}

}  // namespace DRV
}  // namespace NEM
