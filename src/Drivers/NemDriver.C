#include "Drivers/NemDriver.H"

#include "Drivers/DriverJsonTypeTraits.H"

namespace NEM {
namespace DRV {

//---------------------------- Factory of Drivers ----------------------------//
std::unique_ptr<NemDriver> NemDriver::readJSON(const jsoncons::json &inputjson) {
  return inputjson.as<std::unique_ptr<NemDriver>>();
}

DriverOutFile::DriverOutFile(std::string output)
    : outputFile(std::move(output)) {}

DriverInOutFiles::DriverInOutFiles(std::string input, std::string output)
    : inputMeshFile(std::move(input)), outputMeshFile(std::move(output)) {}

}  // namespace DRV
}  // namespace NEM
