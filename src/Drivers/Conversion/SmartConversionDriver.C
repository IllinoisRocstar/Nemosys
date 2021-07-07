#include "Drivers/Conversion/SmartConversionDriver.H"

#include "Mesh/geoMeshFactory.H"

namespace NEM {
namespace DRV {

SmartConversionDriver::SmartConversionDriver(Files files)
    : files_(std::move(files)) {}

SmartConversionDriver::SmartConversionDriver()
    : SmartConversionDriver({{}, {}}) {}

const SmartConversionDriver::Files &SmartConversionDriver::getFiles() const {
  return files_;
}

void SmartConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

void SmartConversionDriver::execute() const {
  vtkSmartPointer<NEM::MSH::geoMeshBase> srcGM =
      NEM::MSH::Read(this->files_.inputMeshFile);
  vtkSmartPointer<NEM::MSH::geoMeshBase> trgGM =
      NEM::MSH::New(this->files_.outputMeshFile);

  trgGM->takeGeoMesh(srcGM);
  trgGM->write(this->files_.outputMeshFile);
}

const SmartConversionDriver::Opts &SmartConversionDriver::getOpts() const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
