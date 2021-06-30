#ifndef NEMOSYS_PROTEUSJSON_H_
#define NEMOSYS_PROTEUSJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/ProteusDriver.H"

// ProteusDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::ProteusDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles, JSONCONS_RDWR),
    (getOpts, setOpts, "Proteus Options", JSONCONS_RDWR),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::ProteusDriver::programType;
     }))

JSONCONS_ALL_MEMBER_NAME_TRAITS(NEM::DRV::ProteusDriver::Files,
                                (fieldFName, "HDF5 Field File"),
                                (meshFName, NEM::DRV::JSON::outMeshFile),
                                (exoMeshFName, "Output Exodus File"))

JSONCONS_N_MEMBER_NAME_TRAITS(NEM::DRV::ProteusDriver::Opts, 1,
                              (edgeSidesetName, "Edge Sideset"),
                              (lowOrder, "Output Low Order"),
                              (bndryConst, "Constraint Boundary"))

#endif  // NEMOSYS_PROTEUSJSON_H_
