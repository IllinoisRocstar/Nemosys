#ifndef NEMOSYS_NEMDRIVER_H_
#define NEMOSYS_NEMDRIVER_H_

#include <jsoncons/json.hpp>

#include "nemosys_export.h"

// TODO: need to harden json input checks etc.
class NEMOSYS_EXPORT NemDriver {
 public:
  NemDriver() = default;
  virtual ~NemDriver() = default;

  static NemDriver *readJSON(const jsoncons::json &inputjson);
};

#endif  // NEMOSYS_NEMDRIVER_H_
