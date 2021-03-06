#ifndef NEMOSYS_MESHQUALITYDRIVER_H_
#define NEMOSYS_MESHQUALITYDRIVER_H_

#include "nemosys_export.h"
#include "NemDriver.H"
#include "meshBase.H"

#include <string>

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT MeshQualityDriver : public NemDriver {
 public:
  MeshQualityDriver() : mesh(nullptr) {}
  MeshQualityDriver(const std::string &_mesh, const std::string &ofname);
  ~MeshQualityDriver() override;

  static MeshQualityDriver *readJSON(const jsoncons::json &inputjson);

 private:
  meshBase *mesh;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_MESHQUALITYDRIVER_H_
