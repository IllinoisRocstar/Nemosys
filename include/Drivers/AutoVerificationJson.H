#ifndef NEMOSYS_AUTOVERIFICATIONJSON_H_
#define NEMOSYS_AUTOVERIFICATIONJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/AutoVerificationDriver.H"

// AutoVerificationDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::AutoVerificationDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getOpts, setOpts, "Verification Options"),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::AutoVerificationDriver::programType;
     }))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::AutoVerificationDriver::Files,
                                (coarseMeshFile, "Coarse Mesh File"),
                                (fineMeshFile, "Fine Mesh File"),
                                (finerMeshFile, "Finer Mesh File"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::AutoVerificationDriver::Opts, 1,
                              (arrayIds, "Array IDs"),
                              (transferType, "Transfer Type"),
                              (targetGCI, "Target GCI"),
                              (numThreads, "Threads"))

#endif  // NEMOSYS_AUTOVERIFICATIONJSON_H_
