#include "Drivers/Conversion/VtkHexToTetConversionDriver.H"

namespace NEM {
namespace DRV {

VtkHexToTetConversionDriver::VtkHexToTetConversionDriver(Files files)
    : files_(std::move(files)) {}

VtkHexToTetConversionDriver::VtkHexToTetConversionDriver()
    : VtkHexToTetConversionDriver({{}, {}}) {}

const VtkHexToTetConversionDriver::Files &
VtkHexToTetConversionDriver::getFiles() const {
  return files_;
}

void VtkHexToTetConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

void VtkHexToTetConversionDriver::execute() const {
  if (this->files_.inputMeshFile.find(".vt") != std::string::npos) {
    std::cout << "Detected file in VTK format" << std::endl;
    std::cout << "Converting to tet Mesh ...." << std::endl;
  } else {
    std::cerr << "Source mesh file is not in VTK format" << std::endl;
  }

  std::ifstream meshStream(this->files_.inputMeshFile);
  if (!meshStream.good()) {
    std::cerr << "Error opening file " << this->files_.inputMeshFile
              << std::endl;
    exit(1);
  }
  // create meshBase object
  std::shared_ptr<meshBase> myMesh =
      meshBase::CreateShared(this->files_.inputMeshFile);

  // Converts hex mesh to tet mesh and writes in VTU file.
  myMesh->convertHexToTetVTK(myMesh->getDataSet());
  myMesh->report();
  myMesh->write(this->files_.outputMeshFile);
}

const VtkHexToTetConversionDriver::Opts &VtkHexToTetConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
