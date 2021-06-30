#ifndef NEMOSYS_MANIPEXOCONVERSIONDRIVER_H_
#define NEMOSYS_MANIPEXOCONVERSIONDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Conversion/ConversionDriver.H"

namespace NEM {
namespace DRV {

/**
 * @brief Manipulates an existing EXODUSII mesh
 */
class NEMOSYS_EXPORT ManipExoConversionDriver : public ConversionDriver {
 public:
  using Files = DriverInOutFiles;

  struct NEMOSYS_EXPORT CombineBlocks {
    std::vector<int> blockIds;
    std::string newName;
  };

  struct NEMOSYS_EXPORT Opts {
    std::vector<CombineBlocks> combineBlocks;
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    static constexpr const char *methodName = "EXO->EXO";
  };

  ManipExoConversionDriver(Files files, Opts opts);

  const Files &getFiles() const;
  void setFiles(Files files);
  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  ManipExoConversionDriver();

  Files files_;
  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_MANIPEXOCONVERSIONDRIVER_H_
