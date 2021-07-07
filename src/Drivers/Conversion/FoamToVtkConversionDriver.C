#include "Drivers/Conversion/FoamToVtkConversionDriver.H"

#include "Mesh/foamMesh.H"
#include "Mesh/vtkMesh.H"

namespace NEM {
namespace DRV {

FoamToVtkConversionDriver::FoamToVtkConversionDriver(Files file)
    : file_(std::move(file)) {}

FoamToVtkConversionDriver::FoamToVtkConversionDriver()
    : FoamToVtkConversionDriver(Files{{}}) {}

const FoamToVtkConversionDriver::Files &FoamToVtkConversionDriver::getFiles()
const {
  return file_;
}

void FoamToVtkConversionDriver::setFiles(Files file) {
  this->file_ = std::move(file);
}

void FoamToVtkConversionDriver::execute() const {
  meshBase *fm = new FOAM::foamMesh();
  fm->read("NULL");
  // TODO: Fix report and write methods for the foamMesh class
  // std::cout << "Variable values is = " << srcmsh << std::endl;
  vtkMesh *vm = new vtkMesh(fm->getDataSet(), this->file_.outputFile);
  vm->report();
  vm->write();
  delete vm;
  delete fm;
}

const FoamToVtkConversionDriver::Opts &FoamToVtkConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
