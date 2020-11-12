#include "NemDriver.H"

#include <iostream>
#include <string>

#include "ConversionDriver.H"
#include "InputGenDriver.H"
#include "MeshGenDriver.H"
#include "MeshQualityDriver.H"
#include "RefineDriver.H"
#include "TransferDriver.H"

#include "NucMeshDriver.H"
#include "PackMeshDriver.H"
#ifdef HAVE_HDF5
#  include "ProteusDriver.H"
#endif
#ifdef HAVE_CGNS
#  include "RocPartCommGenDriver.H"
#endif
#ifdef HAVE_SIMMETRIX
#  include "RemeshDriver.H"
#endif
#ifdef HAVE_TEMPLATE_MESH
#  include "TemplateMeshDriver.H"
#endif

namespace NEM {
namespace DRV {

//---------------------------- Factory of Drivers ----------------------------//
NemDriver *NemDriver::readJSON(const jsoncons::json &inputjson) {
  std::string program_type = inputjson["Program Type"].as<std::string>();
  if (program_type == "Transfer") {
    return TransferDriver::readJSON(inputjson);
  } else if (program_type == "Refinement") {
    return RefineDriver::readJSON(inputjson);
  } else if (program_type == "Mesh Generation") {
    return MeshGenDriver::readJSON(inputjson);
  } else if (program_type == "Mesh Quality") {
    return MeshQualityDriver::readJSON(inputjson);
  } else if (program_type == "Conversion") {
    return ConversionDriver::readJSON(inputjson);
  } else if (program_type == "Input Generation") {
    return InputGenDriver::readJSON(inputjson);
  } else if (program_type == "NucMesh Generation") {
    return NucMeshDriver::readJSON(inputjson);
  } else if (program_type == "Pack Mesh Generation") {
    return PackMeshDriver::readJSON(inputjson);
  } else if (program_type == "Rocstar Remeshing") {
#ifdef HAVE_SIMMETRIX
    return RemeshDriver::readJSON(inputjson);
#else
    std::cerr << "Program Type " << program_type
              << " is not enabled. Build NEMoSys with Simmetrix capabilities."
              << std::endl;
    exit(1);
#endif  // HAVE_SIMMETRIX
    /*
  } else if (program_type == "Post Rocstar Remeshing") {
    return RocRestartDriver::readJSON(inputjson);
  } else if (program_type == "Rocstar Communication Generation") {
    return RocPrepDriver::readJSON(inputjson);
    */
  } else if (program_type == "Rocstar Communication Generation") {
#ifdef HAVE_CGNS
    return RocPartCommGenDriver::readJSON(inputjson);
#else
    std::cerr << "Program Type " << program_type
              << " is not enabled. Build NEMoSys with CGNS capabilities."
              << std::endl;
    exit(1);
#endif  // HAVE_CGNS
  } else if (program_type == "Template Mesh Generation") {
#ifdef HAVE_TEMPLATE_MESH
    return TemplateMeshDriver::readJSON(inputjson);
#else
    std::cerr
        << "Program Type " << program_type
        << " is not enabled. Build NEMoSys with template mesh capabilities."
        << std::endl;
    exit(1);
#endif  // HAVE_TEMPLATE_MESH
  } else if (program_type == "Proteus") {
#ifdef HAVE_HDF5
    return ProteusDriver::readJSON(inputjson);
#else
    std::cerr << "Program Type " << program_type
              << " is not enabled. Build NEMoSys with HDF5 capabilities."
              << std::endl;
    exit(1);
#endif  // HAVE_HDF5
  } else {
    std::cerr << "Program Type " << program_type
              << " is not supported by NEMoSys." << std::endl;
    exit(1);
  }
}

}  // namespace DRV
}  // namespace NEM
