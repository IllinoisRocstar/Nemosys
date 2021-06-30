#ifndef NEMOSYS_FOAMTOVTKCONVERSIONDRIVER_H_
#define NEMOSYS_FOAMTOVTKCONVERSIONDRIVER_H_

#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT FoamToVtkConversionDriver : public ConversionDriver {
 public:
  using Files = DriverOutFile;

  explicit FoamToVtkConversionDriver(Files file);

  const Files &getFiles() const;
  void setFiles(Files file);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  FoamToVtkConversionDriver();

  Files file_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "FOAM->VTK";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_FOAMTOVTKCONVERSIONDRIVER_H_
