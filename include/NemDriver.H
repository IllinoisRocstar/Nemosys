#ifndef NEMDRIVERS_H
#define NEMDRIVERS_H

#include <meshBase.H>
#include <jsoncons/json.hpp>

using jsoncons::json;

class NemDriver
{
  public:
    NemDriver(){};
    virtual ~NemDriver(){}
    static NemDriver* readJSON(json inputjson);
};

#endif
