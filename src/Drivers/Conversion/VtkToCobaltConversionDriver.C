#include "Drivers/Conversion/VtkToCobaltConversionDriver.H"

#include "Mesh/cobalt.H"

namespace NEM {
namespace DRV {

VtkToCobaltConversionDriver::Files::Files(std::string input,
                                          std::string outputCgr,
                                          std::string outputCgi)
    : inputMeshFile(std::move(input)),
      outputCgrFile(std::move(outputCgr)),
      outputCgiFile(std::move(outputCgi)) {}

VtkToCobaltConversionDriver::VtkToCobaltConversionDriver(Files files)
    : files_(std::move(files)) {}

VtkToCobaltConversionDriver::VtkToCobaltConversionDriver()
    : VtkToCobaltConversionDriver({{}, {}, {}}) {}

const VtkToCobaltConversionDriver::Files &
VtkToCobaltConversionDriver::getFiles() const {
  return files_;
}

void VtkToCobaltConversionDriver::setFiles(Files files) {
  this->files_ = std::move(files);
}

void VtkToCobaltConversionDriver::execute() const {
  if (this->files_.inputMeshFile.find(".vt") != std::string::npos) {
    std::cout << "Detected file in VTK format" << std::endl;
    std::cout << "Converting to COBALT ...." << std::endl;
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

  // create Cobalt object from meshBase
  COBALT::cobalt *cm = new COBALT::cobalt(myMesh, this->files_.inputMeshFile,
                                          this->files_.outputCgrFile,
                                          this->files_.outputCgiFile);
  // write to file
  cm->write();
}

const VtkToCobaltConversionDriver::Opts &VtkToCobaltConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
