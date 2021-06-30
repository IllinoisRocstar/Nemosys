#ifndef NEMOSYS_UNIFORMREFINEDRIVER_H_
#define NEMOSYS_UNIFORMREFINEDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Refine/RefineDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT UniformRefineDriver : public RefineDriver {
 public:
  struct NEMOSYS_EXPORT Opts {
   public:
    Opts(bool transferData, double edgeScale);
    bool transferData{};
    double edgeScale{};
    static constexpr const char *method = "uniform";
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };
  UniformRefineDriver(Files files, Opts opts);

  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  UniformRefineDriver();

  Opts opts_;
};
}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_UNIFORMREFINEDRIVER_H_
