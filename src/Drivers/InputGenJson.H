#ifndef NEMOSYS_INPUTGENJSON_H_
#define NEMOSYS_INPUTGENJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/InputGenDriver.H"

// InputGenDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::InputGenDriver, NEM::DRV::NemDriver, 2,
    (getService, setService, "Service"),
    (getOpts, setOpts, "Input Generation Options"),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::InputGenDriver::programType;
     }))

#endif  // NEMOSYS_INPUTGENJSON_H_
