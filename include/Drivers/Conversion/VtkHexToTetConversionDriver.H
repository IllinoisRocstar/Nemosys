#ifndef NEMOSYS_VTKHEXTOTETCONVERSIONDRIVER_H_
#define NEMOSYS_VTKHEXTOTETCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT VtkHexToTetConversionDriver : public ConversionDriver {
 public:
  using Files = DriverInOutFiles;

  explicit VtkHexToTetConversionDriver(Files files);

  const Files &getFiles() const;
  void setFiles(Files files);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  VtkHexToTetConversionDriver();

  Files files_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "VTK-HEX->VTK-TET";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_VTKHEXTOTETCONVERSIONDRIVER_H_
