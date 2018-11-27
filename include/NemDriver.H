#ifndef NEMDRIVERS_H
#define NEMDRIVERS_H

#include <meshBase.H>
#include <jsoncons/json.hpp>

using jsoncons::json;

// TODO: need to harden json input checks etc.

class NemDriver
{
  public:
    NemDriver(){}
    virtual ~NemDriver(){}
    static NemDriver* readJSON(json inputjson);
};

#endif
