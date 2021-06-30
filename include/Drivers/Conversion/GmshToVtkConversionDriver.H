#ifndef NEMOSYS_GMSHTOVTKCONVERSIONDRIVER_H_
#define NEMOSYS_GMSHTOVTKCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT GmshToVtkConversionDriver : public ConversionDriver {
 public:
  using Files = DriverInOutFiles;

  explicit GmshToVtkConversionDriver(Files files);

  const Files &getFiles() const;
  void setFiles(Files files);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  GmshToVtkConversionDriver();

  Files files_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "GMSH->VTK";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_GMSHTOVTKCONVERSIONDRIVER_H_
