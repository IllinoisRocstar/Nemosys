#include "Drivers/PackMesh/PackMeshDriver.H"

namespace NEM {
namespace DRV {

jsoncons::string_view PackMeshDriver::getProgramType() const {
  return programType;
}

}  // namespace DRV
}  // namespace NEM
