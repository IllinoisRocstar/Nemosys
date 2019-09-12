#ifndef INPUTGENDRIVER_H
#define INPUTGENDRIVER_H

#include "nemosys_export.h"
#include "NemDriver.H"

#include <string>

class NEMOSYS_EXPORT InputGenDriver : public NemDriver
{
  public:
    InputGenDriver(const std::string &srvName,
                   const jsoncons::json &inputjson);
    ~InputGenDriver() override;

    static InputGenDriver *readJSON(const jsoncons::json &inputjson);
    static InputGenDriver *readJSON(const std::string &ifname);
};

#endif
