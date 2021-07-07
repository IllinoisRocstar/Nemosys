#include "Drivers/Conversion/FoamToMshConversionDriver.H"

#include "Mesh/foamMesh.H"
#include "Mesh/gmshMesh.H"

namespace NEM {
namespace DRV {

FoamToMshConversionDriver::FoamToMshConversionDriver(Files file)
    : file_(std::move(file)) {}

FoamToMshConversionDriver::FoamToMshConversionDriver()
    : FoamToMshConversionDriver(DriverOutFile{{}}) {}

const FoamToMshConversionDriver::Files &FoamToMshConversionDriver::getFiles()
    const {
  return file_;
}

void FoamToMshConversionDriver::setFiles(Files file) {
  this->file_ = std::move(file);
}

void FoamToMshConversionDriver::execute() const {
  meshBase *fm = new FOAM::foamMesh();
  fm->read("NULL");
  // TODO: Fix report and write methods for the foamMesh class
  // fm->setFileName(ofname);
  // fm->report();
  // fm->writeMSH();
  auto *gm = new gmshMesh(fm);
  gm->write(this->file_.outputFile);
  delete fm;
}

const FoamToMshConversionDriver::Opts &FoamToMshConversionDriver::getOpts()
    const {
  static constexpr Opts opts{};
  return opts;
}

}  // namespace DRV
}  // namespace NEM
