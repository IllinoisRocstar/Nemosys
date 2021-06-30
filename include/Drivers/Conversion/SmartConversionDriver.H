#ifndef NEMOSYS_SMARTCONVERSIONDRIVER_H_
#define NEMOSYS_SMARTCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT SmartConversionDriver : public ConversionDriver {
 public:
  using Files = DriverInOutFiles;

  explicit SmartConversionDriver(Files files);

  const Files &getFiles() const;
  void setFiles(Files files);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  SmartConversionDriver();

  Files files_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "SMART";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_SMARTCONVERSIONDRIVER_H_
