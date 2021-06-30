#include "Drivers/Conversion/VtkToPntConversionDriver.H"

namespace NEM {
namespace DRV {

VtkToPntConversionDriver::Opts::Opts(int dim, PNTMesh::BlockMap blockMap)
    : dim(dim), elemBlockMap(std::move(blockMap)) {}

VtkToPntConversionDriver::VtkToPntConversionDriver(Files files, Opts opts)
    : files_(std::move(files)), opts_(std::move(opts)) {}

VtkToPntConversionDriver::VtkToPntConversionDriver()
    : VtkToPntConversionDriver({{}, {}}, {{}, {}}) {}

const VtkToPntConversionDriver::Files &VtkToPntConversionDriver::getFiles()
    const {
  return files_;
}

void VtkToPntConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

const VtkToPntConversionDriver::Opts &VtkToPntConversionDriver::getOpts()
    const {
  return opts_;
}

void VtkToPntConversionDriver::setOpts(Opts opts) {
  this->opts_ = std::move(opts);
}

void VtkToPntConversionDriver::execute() const {
  auto source = meshBase::Create(files_.inputMeshFile);
  std::cout << "Number of Blocks : " << opts_.elemBlockMap.size() << std::endl;
  auto *pm = new PNTMesh::pntMesh(source, opts_.dim, opts_.elemBlockMap.size(),
                                  opts_.elemBlockMap);
  pm->write(files_.outputMeshFile);
  delete pm;
}

}  // namespace DRV
}  // namespace NEM
