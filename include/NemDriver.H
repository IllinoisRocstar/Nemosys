#ifndef NEMDRIVERS_H
#define NEMDRIVERS_H

#include "nemosys_export.h"
#include "meshBase.H"

#include <jsoncons/json.hpp>

// TODO: need to harden json input checks etc.
class NEMOSYS_EXPORT NemDriver {
  public:
    NemDriver() = default;
    virtual ~NemDriver() = default;

    static NemDriver *readJSON(const jsoncons::json &inputjson);
};

#endif
