#include "Drivers/MeshGen/MeshGenDriver.H"

namespace NEM {
namespace DRV {

MeshGenDriver::MeshGenFiles::MeshGenFiles(std::string geometry,
                                          std::string output)
    : inputGeoFile(std::move(geometry)), outputMeshFile(std::move(output)) {}

jsoncons::string_view MeshGenDriver::getProgramType() const {
  return programType;
}

}  // namespace DRV
}  // namespace NEM
