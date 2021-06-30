#ifndef NEMOSYS_VTKTOPNTCONVERSIONDRIVER_H_
#define NEMOSYS_VTKTOPNTCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

#include "pntMesh.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT VtkToPntConversionDriver : public ConversionDriver {
 public:
  using Files = DriverInOutFiles;

  struct NEMOSYS_EXPORT Opts {
   public:
    Opts(int dim, PNTMesh::BlockMap blockMap);
    int dim{};
    PNTMesh::BlockMap elemBlockMap{};
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *methodName = "VTK->PNT";
  };

  VtkToPntConversionDriver(Files files, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files files);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  VtkToPntConversionDriver();

  Files files_;
  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_VTKTOPNTCONVERSIONDRIVER_H_
