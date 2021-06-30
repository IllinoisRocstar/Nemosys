#include "Drivers/Conversion/ManipExoConversionDriver.H"

#include "exoMesh.H"

namespace NEM {
namespace DRV {

ManipExoConversionDriver::ManipExoConversionDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

ManipExoConversionDriver::ManipExoConversionDriver()
    : ManipExoConversionDriver(Files{{}, {}}, Opts{{}}) {}

const ManipExoConversionDriver::Files &ManipExoConversionDriver::getFiles()
    const {
  return files_;
}

void ManipExoConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const ManipExoConversionDriver::Opts &ManipExoConversionDriver::getOpts()
    const {
  return opts_;
}

void ManipExoConversionDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void ManipExoConversionDriver::execute() const {
  // read in exo mesh from file name
  NEM::MSH::EXOMesh::exoMesh em{this->files_.inputMeshFile};
  em.read(this->files_.inputMeshFile);

  // Combining Blocks
  for (const auto &cmb : this->opts_.combineBlocks) {
    em.combineElmBlks(cmb.blockIds, cmb.newName);
  }

  // Set the output file name
  em.setFileName(this->files_.outputMeshFile);
  em.write();
  em.report();
}

}  // namespace DRV
}  // namespace NEM
