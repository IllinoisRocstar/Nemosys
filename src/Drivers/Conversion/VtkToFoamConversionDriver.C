#include "Drivers/Conversion/VtkToFoamConversionDriver.H"

#include "Mesh/foamMesh.H"

namespace NEM {
namespace DRV {

VtkToFoamConversionDriver::Files::Files(std::string input)
    : inputMeshFile(std::move(input)) {}

VtkToFoamConversionDriver::VtkToFoamConversionDriver(Files file)
    : file_(std::move(file)) {}

VtkToFoamConversionDriver::VtkToFoamConversionDriver()
    : VtkToFoamConversionDriver(Files{{}}) {}

const VtkToFoamConversionDriver::Files &VtkToFoamConversionDriver::getFiles()
    const {
  return file_;
}

void VtkToFoamConversionDriver::setFiles(Files file) {
  this->file_ = std::move(file);
}

void VtkToFoamConversionDriver::execute() const {
  if (this->file_.inputMeshFile.find(".vt") != std::string::npos) {
    std::cout << "Detected file in VTK format" << std::endl;
    std::cout << "Converting to FOAM Mesh ...." << std::endl;
  } else {
    std::cerr << "Source mesh file is not in VTK format" << std::endl;
  }

  std::ifstream meshStream(this->file_.inputMeshFile);
  if (!meshStream.good()) {
    std::cerr << "Error opening file " << this->file_.inputMeshFile
              << std::endl;
    exit(1);
  }
  // create meshBase object
  std::shared_ptr<meshBase> myMesh =
      meshBase::CreateShared(this->file_.inputMeshFile);

  // create foamMesh object
  FOAM::foamMesh *fm = new FOAM::foamMesh(myMesh);

  // Write polyMesh
  // fm->report();
  fm->write({});
}

const VtkToFoamConversionDriver::Opts &VtkToFoamConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
