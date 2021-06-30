#ifndef NEMOSYS_VTKTOCOBALTCONVERSIONDRIVER_H_
#define NEMOSYS_VTKTOCOBALTCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT VtkToCobaltConversionDriver : public ConversionDriver {
 public:
  struct NEMOSYS_EXPORT Files {
   public:
    Files(std::string input, std::string outputCgr, std::string outputCgi);
    std::string inputMeshFile{};
    std::string outputCgrFile{};
    std::string outputCgiFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Files() = default;
  };

  explicit VtkToCobaltConversionDriver(Files files);

  const Files &getFiles() const;
  void setFiles(Files files);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  VtkToCobaltConversionDriver();

  Files files_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "VTK->COBALT";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_VTKTOCOBALTCONVERSIONDRIVER_H_
