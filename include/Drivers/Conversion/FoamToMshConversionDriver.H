#ifndef NEMOSYS_FOAMTOMSHCONVERSIONDRIVER_H_
#define NEMOSYS_FOAMTOMSHCONVERSIONDRIVER_H_

#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT FoamToMshConversionDriver : public ConversionDriver {
 public:
  using Files = DriverOutFile;

  explicit FoamToMshConversionDriver(Files file);

  const Files &getFiles() const;
  void setFiles(Files file);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  FoamToMshConversionDriver();

  Files file_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "FOAM->MSH";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_FOAMTOMSHCONVERSIONDRIVER_H_
