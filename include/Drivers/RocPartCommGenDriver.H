#ifndef NEMOSYS_ROCPARTCOMMGENDRIVER_H_
#define NEMOSYS_ROCPARTCOMMGENDRIVER_H_

#include <unordered_set>

#include <vtkEdgeTable.h>
#include <vtkPolyData.h>

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"
#include "Mesh/meshBase.H"

#include <cgnslib.h>

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT RocPartCommGenDriver : public NemDriver {
 public:
  RocPartCommGenDriver(std::string volName, std::string surfName,
                       int numPartitions);

  std::string volName{};
  std::string surfName{};
  int numPartitions{};

  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  RocPartCommGenDriver() = default;

  static constexpr const char *programType = "Rocstar Communication Generation";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_ROCPARTCOMMGENDRIVER_H_
