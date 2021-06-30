#ifndef NEMOSYS_TEMPLATEMESHJSON_H_
#define NEMOSYS_TEMPLATEMESHJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/TemplateMeshDriver.H"

// TemplateMeshDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::TemplateMeshDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, "Template Options", JSONCONS_RDWR),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::TemplateMeshDriver::programType;
     }))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::TemplateMeshDriver::Opts, 1,
                              (templateName, "Template Name"),
                              (templateParams, "Params"))

#endif  // NEMOSYS_TEMPLATEMESHJSON_H_
