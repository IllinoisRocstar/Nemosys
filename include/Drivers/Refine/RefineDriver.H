#ifndef NEMOSYS_REFINEDRIVER_H_
#define NEMOSYS_REFINEDRIVER_H_

#include <string>

#include "nemosys_export.h"
#include "Drivers/NemDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT RefineDriver : public NemDriver {
 public:
  using Files = DriverInOutFiles;
  explicit RefineDriver(Files files);

  virtual const Files &getFiles() const;
  virtual void setFiles(Files files);

  JSONCONS_TYPE_TRAITS_FRIEND

 protected:
  Files files_;

 private:
  static constexpr const char *programType = "Refinement";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_REFINEDRIVER_H_
