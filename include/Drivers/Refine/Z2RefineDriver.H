#ifndef NEMOSYS_Z2REFINEDRIVER_H_
#define NEMOSYS_Z2REFINEDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Refine/RefineDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT Z2RefineDriver : public RefineDriver {
 public:
  struct NEMOSYS_EXPORT Opts {
   public:
    Opts(bool transferData, std::string arrayName, int order);
    bool transferData{};
    std::string arrayName{};
    int order{};
    static constexpr const char *method = "Z2 Error Estimator";
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
  };
  Z2RefineDriver(Files files, Opts opts);

  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  Z2RefineDriver();

  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_Z2REFINEDRIVER_H_
