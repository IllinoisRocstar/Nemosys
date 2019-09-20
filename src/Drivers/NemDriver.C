#include "NemDriver.H"

#include <string>
#include <iostream>

#include "TransferDriver.H"
#include "RefineDriver.H"
#include "MeshQualityDriver.H"
#include "MeshGenDriver.H"
#include "ConversionDriver.H"
#include "InputGenDriver.H"

#include "NucMeshDriver.H"
#include "PackMeshDriver.H"
#ifdef HAVE_CGNS
#  include "RocPartCommGenDriver.H"
#endif
#ifdef HAVE_SIMMETRIX
#  include "RemeshDriver.H"
#endif



//------------------------------ Factory of Drivers ----------------------------------------//
NemDriver *NemDriver::readJSON(const jsoncons::json &inputjson)
{
  std::string program_type = inputjson["Program Type"].as<std::string>();
  if (program_type == "Transfer")
  {
    return TransferDriver::readJSON(inputjson);
  }
  else if (program_type == "Refinement")
  {
    return RefineDriver::readJSON(inputjson);
  }
  else if (program_type == "Mesh Generation")
  {
    return MeshGenDriver::readJSON(inputjson);
  }
  else if (program_type == "Mesh Quality")
  {
    return MeshQualityDriver::readJSON(inputjson);
  }
  else if (program_type == "Conversion")
  {
    return ConversionDriver::readJSON(inputjson);
  }
  else if (program_type == "Input Generation")
  {
    return InputGenDriver::readJSON(inputjson);
  }
  else if (!program_type.compare("NucMesh Generation"))
  {
  	return NucMeshDriver::readJSON(inputjson);
  }
  else if (program_type == "Pack Mesh Generation")
  {
    #ifdef HAVE_CFMSH
    return PackMeshDriver::readJSON(inputjson);
    #else
    std::cerr << "Build NEMoSys with CfMesh" << std::endl;
    exit(1);
    #endif
  }
  else if (program_type == "Rocstar Remeshing")
  {
#ifdef HAVE_SIMMETRIX
    return RemeshDriver::readJSON(inputjson);
#else
    std::cerr << "Program Type " << program_type
              << " is not enabled. Build NEMoSys with Simmetrix capabilities."
              << std::endl;
    exit(1);
#endif // HAVE_SIMMETRIX
  }
  //else if (program_type == "Post Rocstar Remeshing")
  //{
  //  return RocRestartDriver::readJSON(inputjson);
  //}
  //else if (program_type == "Rocstar Communication Generation")
  //{
  //  return RocPrepDriver::readJSON(inputjson);
  //}
  else if (program_type == "Rocstar Communication Generation")
  {
#ifdef HAVE_CGNS
    return RocPartCommGenDriver::readJSON(inputjson);
#else
    std::cerr << "Program Type " << program_type
              << " is not enabled. Build NEMoSys with CGNS capabilities."
              << std::endl;
    exit(1);
#endif  // HAVE_CGNS
  }
  else
  {
    std::cerr << "Program Type " << program_type
              << " is not supported by NEMoSys." << std::endl;
    exit(1);
  }
}
