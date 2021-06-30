#ifndef INPUTGENDRIVER_H
#define INPUTGENDRIVER_H

#include "nemosys_export.h"
#include "NemDriver.H"

namespace NEM {
namespace DRV {

class NEMOSYS_EXPORT InputGenDriver : public NemDriver {
 public:
  InputGenDriver(std::string service, jsoncons::json opts);

  const std::string &getService() const;
  void setService(std::string service);
  const jsoncons::json &getOpts() const;
  void setOpts(jsoncons::json opts);
  void execute() const override;

  JSONCONS_TYPE_TRAITS_FRIEND

 private:
  InputGenDriver() = default;

  static constexpr const char *programType = "Input Generation";

  std::string service_{};
  jsoncons::json opts_{};

  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif
