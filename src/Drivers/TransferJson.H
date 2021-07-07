#ifndef NEMOSYS_TRANSFERJSON_H_
#define NEMOSYS_TRANSFERJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/TransferDriver.H"

// TransferDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::TransferDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, "Transfer Options", JSONCONS_RDWR),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::TransferDriver::programType;
     }))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::TransferDriver::Files,
                                (sourceMeshFile, "Source Mesh File"),
                                (targetMeshFile, "Target Mesh File"),
                                (outputMeshFile, NEM::DRV::JSON::outMeshFile))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::TransferDriver::Opts, 2,
                              (method, "Method"),
                              (checkQuality, "Check Transfer Quality"),
                              (arrayNames, "Array Names"))

#endif  // NEMOSYS_TRANSFERJSON_H_
