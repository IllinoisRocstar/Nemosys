#include "Drivers/Conversion/GmshToVtkConversionDriver.H"

namespace NEM {
namespace DRV {

GmshToVtkConversionDriver::GmshToVtkConversionDriver(Files files)
    : files_(std::move(files)) {}

GmshToVtkConversionDriver::GmshToVtkConversionDriver()
    : GmshToVtkConversionDriver({{}, {}}) {}

const GmshToVtkConversionDriver::Files &GmshToVtkConversionDriver::getFiles()
    const {
  return files_;
}

void GmshToVtkConversionDriver::setFiles(Files files) {
  files_ = std::move(files);
}

void GmshToVtkConversionDriver::execute() const {
  if (this->files_.inputMeshFile.find(".msh") != std::string::npos) {
    std::cout << "Detected file in GMSH format" << std::endl;
    std::cout << "Converting to VTK ...." << std::endl;
  } else {
    std::cerr << "Source mesh file is not in GMSH format" << std::endl;
  }
  meshBase *mb = meshBase::exportGmshToVtk(this->files_.inputMeshFile);
  mb->report();
  mb->write(this->files_.outputMeshFile);
}

const GmshToVtkConversionDriver::Opts &GmshToVtkConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
