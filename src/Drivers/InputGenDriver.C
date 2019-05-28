// misc headers
#include <iostream>
#include <algorithm>
#include <map>

// Nemosys headers
#include "AuxiliaryFunctions.H"
#include "InputGenDriver.H"
//#include "meshSrch.H"


#ifdef HAVE_EPIC
#include "ep16Prep.H"
#endif

//----------------------- Conversion Driver -----------------------------------------//
InputGenDriver::InputGenDriver(std::string srvName, json inputjson)
{
  std::cout << "InputGenDriver created" << std::endl;
  srvName = toLower(srvName);
  Timer T;
  T.start();
  if (!srvName.compare("epic_2016"))
  {
#ifdef HAVE_EPIC
    ep16Prep::readJSON(inputjson);
#else
    std::cerr << "Compile the code with EPIC module enabled.\n";
    exit(0);
#endif  
  }
  T.stop();
}

InputGenDriver::~InputGenDriver()
{
  std::cout << "InputGenDriver destroyed" << std::endl;
}

InputGenDriver* InputGenDriver::readJSON(json inputjson)
{
  std::string srvName = inputjson["Service"].as<std::string>();
  InputGenDriver* inpGenDrv;
  inpGenDrv = new InputGenDriver(srvName, inputjson);
  return inpGenDrv;
}

InputGenDriver* InputGenDriver::readJSON(std::string ifname)
{
  std::cout << "Reading JSON file\n";
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || find_ext(ifname) != ".json")
  {
    std::cout << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (find_ext(ifname) != ".json")
  {
    std::cout << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  json inputjson;
  inputStream >> inputjson;
  
  // checking if array
  if (inputjson.is_array())
  {
    std::cout << "Warning: Input is an array. Only first element will be processed\n";
    return InputGenDriver::readJSON(inputjson[0]);
  } 
  else
  {
    return InputGenDriver::readJSON(inputjson);
  }
    
}

