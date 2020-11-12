#ifndef NEMOSYS_NEMDRIVER_H_
#define NEMOSYS_NEMDRIVER_H_

#include <jsoncons/json.hpp>

#include "nemosys_export.h"

namespace NEM {
namespace DRV {

// TODO: need to harden json input checks etc.
/**
 * @class NemDriver
 * @brief base class for drivers
 */
class NEMOSYS_EXPORT NemDriver {
 public:
  NemDriver() = default;
  virtual ~NemDriver() = default;

  /**
   * Factory method for all drivers. Calls derived class's readJSON method.
   * @param inputjson json input with "Program Type" key
   * @return pointer to base NemDriver class
   */
  static NemDriver *readJSON(const jsoncons::json &inputjson);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_NEMDRIVER_H_
