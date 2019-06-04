#include <NemDriver.H>
#include <TransferDriver.H>
#include <MeshQualityDriver.H>
#include <MeshGenDriver.H>
#include <RefineDriver.H>
#include <ConversionDriver.H>
#include "InputGenDriver.H"
#include <RemeshDriver.H>
#include <RocPartCommGenDriver.H>

//------------------------------ Factory of Drivers ----------------------------------------//
NemDriver* NemDriver::readJSON(json inputjson)
{
  std::string program_type = inputjson["Program Type"].as<std::string>();
  if (!program_type.compare("Transfer"))
  {
    return TransferDriver::readJSON(inputjson); 
  }
  else if (!program_type.compare("Refinement"))
  {
    return RefineDriver::readJSON(inputjson);
  }
  else if (!program_type.compare("Mesh Generation"))
  {
    return MeshGenDriver::readJSON(inputjson);
  }
  else if (!program_type.compare("Mesh Quality"))
  {
    return MeshQualityDriver::readJSON(inputjson);
  }
  else if (!program_type.compare("Conversion"))
  {
    return ConversionDriver::readJSON(inputjson);
  }
  else if (!program_type.compare("Input Generation"))
  {
    return InputGenDriver::readJSON(inputjson);
  }
  else if (!program_type.compare("Rocstar Remeshing"))
  {
#ifdef HAVE_SIMMETRIX
    return RemeshDriver::readJSON(inputjson);
#else
    std::cerr << "Program Type " << program_type
              << " is not enabled. Build NEMoSys with Simmetrix capabilities." << std::endl;
    exit(1);
#endif // HAVE_SIMMETRIX
  }
  //else if (!program_type.compare("Post Rocstar Remeshing"))
  //{
  //  return RocRestartDriver::readJSON(inputjson);
  //}
  //else if (!program_type.compare("Rocstar Communication Generation"))
  //{
  //  return RocPrepDriver::readJSON(inputjson);
  //}
  else if (!program_type.compare("Rocstar Communication Generation"))
  {
    return RocPartCommGenDriver::readJSON(inputjson);
  }
  else
  {
    std::cerr << "Program Type " << program_type
              << " is not supported by NEMoSys." << std::endl;
    exit(1);
  }
}

