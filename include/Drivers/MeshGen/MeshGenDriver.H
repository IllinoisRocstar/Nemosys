#ifndef NEMOSYS_MESHGENDRIVER_H_
#define NEMOSYS_MESHGENDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT MeshGenDriver : public NemDriver {
 public:
  struct NEMOSYS_EXPORT MeshGenFiles {
   public:
    MeshGenFiles(std::string geometry, std::string output);

    std::string inputGeoFile{};
    std::string outputMeshFile{};

    JSONCONS_TYPE_TRAITS_FRIEND

   private:
    MeshGenFiles() = default;
  };

  JSONCONS_TYPE_TRAITS_FRIEND

 protected:
  MeshGenDriver() = default;

 private:
  static constexpr const char *programType = "Mesh Generation";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_MESHGENDRIVER_H_
