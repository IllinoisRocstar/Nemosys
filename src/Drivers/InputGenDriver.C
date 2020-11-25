// Nemosys headers
#include "AuxiliaryFunctions.H"
#include "InputGenDriver.H"


#include <fstream>
#include <iostream>

namespace NEM {
namespace DRV {

//----------------------- Conversion Driver
//-----------------------------------------//
InputGenDriver::InputGenDriver(const std::string &_srvName,
                               const jsoncons::json &inputjson) {
  std::cout << "InputGenDriver created" << std::endl;
  std::string srvName = _srvName;
  nemAux::toLower(srvName);
  nemAux::Timer T;
  T.start();
  T.stop();
}

InputGenDriver::~InputGenDriver() {
  std::cout << "InputGenDriver destroyed" << std::endl;
}

InputGenDriver *InputGenDriver::readJSON(const jsoncons::json &inputjson) {
  std::string srvName = inputjson["Service"].as<std::string>();
  InputGenDriver *inpGenDrv;
  inpGenDrv = new InputGenDriver(srvName, inputjson);
  return inpGenDrv;
}

InputGenDriver *InputGenDriver::readJSON(const std::string &ifname) {
  std::cout << "Reading JSON file" << std::endl;
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  // checking if array
  if (inputjson.is_array()) {
    std::cerr
        << "Warning: Input is an array. Only first element will be processed"
        << std::endl;
    return InputGenDriver::readJSON(inputjson[0]);
  } else {
    return InputGenDriver::readJSON(inputjson);
  }
}

}  // namespace DRV
}  // namespace NEM
