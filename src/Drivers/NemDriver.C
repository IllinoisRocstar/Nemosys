#include <NemDriver.H>
#include <TransferDriver.H>
#include <MeshQualityDriver.H>
#include <MeshGenDriver.H>
#include <RefineDriver.H>
#include <ConversionDriver.H>
#include <RemeshDriver.H>

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
  else if (!program_type.compare("Rocstar Remeshing"))
  {
    return RemeshDriver::readJSON(inputjson);
  }
  else
  {
    std::cout << "Program Type " << program_type 
              << " is not supported by Nemosys" << std::endl;
    exit(1);
  }
  
}
