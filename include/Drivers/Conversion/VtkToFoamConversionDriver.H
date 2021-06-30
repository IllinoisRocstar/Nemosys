#ifndef NEMOSYS_VTKTOFOAMCONVERSIONDRIVER_H_
#define NEMOSYS_VTKTOFOAMCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT VtkToFoamConversionDriver : public ConversionDriver {
 public:
  struct NEMOSYS_EXPORT Files {
   public:
    explicit Files(std::string input);
    std::string inputMeshFile{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Files() = default;
  };

  explicit VtkToFoamConversionDriver(Files file);

  const Files &getFiles() const;
  void setFiles(Files file);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  VtkToFoamConversionDriver();

  Files file_;

 public:
  struct NEMOSYS_NO_EXPORT Opts {
    static constexpr const char *methodName = "VTK->FOAM";
  };

 private:
  const Opts &getOpts() const;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_VTKTOFOAMCONVERSIONDRIVER_H_
