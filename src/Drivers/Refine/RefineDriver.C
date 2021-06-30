#include "Drivers/Refine/RefineDriver.H"

namespace NEM {
namespace DRV {

RefineDriver::RefineDriver(Files files) : files_(std::move(files)) {}

jsoncons::string_view RefineDriver::getProgramType() const {
  return programType;
}

const RefineDriver::Files &RefineDriver::getFiles() const {
  return files_;
}

void RefineDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

}  // namespace DRV
}  // namespace NEM
