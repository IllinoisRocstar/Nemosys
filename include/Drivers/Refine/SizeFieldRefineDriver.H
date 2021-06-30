#ifndef NEMOSYS_VALUEREFINEDRIVER_H_
#define NEMOSYS_VALUEREFINEDRIVER_H_

#include "nemosys_export.h"
#include "Drivers/Refine/RefineDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT SizeFieldRefineDriver : public RefineDriver {
 public:
  struct NEMOSYS_EXPORT Opts {
   public:
    enum class Method { VALUE, GRADIENT };
    Opts(Method method, std::string arrayName, double stdDevMult,
         bool maxIsMin, bool transferData);
    Method method{};
    std::string arrayName{};
    double stdDevMult{};
    bool maxIsMin{};
    bool transferData{};
    double sizeFactor{1.};
    std::string getMethodStr() const;
    JSONCONS_TYPE_TRAITS_FRIEND
   private:
    Opts() = default;
    static constexpr const char *valStr = "value";
    static constexpr const char *gradStr = "gradient";
  };

  SizeFieldRefineDriver(Files files, Opts opts);

  const Opts &getOpts() const;
  void setOpts(Opts opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  SizeFieldRefineDriver();

  Opts opts_;
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_VALUEREFINEDRIVER_H_
