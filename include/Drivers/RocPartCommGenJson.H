#ifndef NEMOSYS_ROCPARTCOMMGENJSON_H_
#define NEMOSYS_ROCPARTCOMMGENJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/RocPartCommGenDriver.H"

// RocPartCommGenDriver
NEM_JSON_N_MEMBER_NAME_TRAITS_FINAL(
    NEM::DRV::RocPartCommGenDriver, NEM::DRV::NemDriver, 3,
    (volName, "Remeshed Volume"), (surfName, "Stitched Surface"),
    (numPartitions, "Number Of Partitions"),
    (programType, NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::RocPartCommGenDriver::programType;
     }))

#endif  // NEMOSYS_ROCPARTCOMMGENJSON_H_
