#ifndef NEMOSYS_NUCMESHJSON_H_
#define NEMOSYS_NUCMESHJSON_H_

#include <jsoncons/json.hpp>
#include "Drivers/NemJsonMacros.H"

#include "Drivers/NucMeshDriver.H"

// NucMeshDriver
NEM_JSON_N_GETTER_SETTER_NAME_TRAITS_FINAL(
    NEM::DRV::NucMeshDriver, NEM::DRV::NemDriver, 2,
    (getFiles, setFiles, NEM::DRV::JSON::meshFiles),
    (getGeometryAndMesh, setGeometryAndMesh, "Geometry and Mesh"),
    (getProgramType, , NEM::DRV::JSON::programType, NEM_JSON_RDONLY_OVERRIDE,
     [](const jsoncons::string_view &x) {
       return x == NEM::DRV::NucMeshDriver::programType;
     }))

#endif  // NEMOSYS_NUCMESHJSON_H_
