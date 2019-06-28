#include "TransferDriver.H"
#include "AuxiliaryFunctions.H"

#include <iostream>
#include <string>

//----------------------- Transfer Driver ------------------------------------//
TransferDriver::TransferDriver(const std::string &srcmsh,
                               const std::string &trgmsh,
                               const std::string &method,
                               const std::string &ofname,
                               bool checkQuality)
{
  source = meshBase::Create(srcmsh);
  target = meshBase::Create(trgmsh);
  std::cout << "TransferDriver created" << std::endl;
  nemAux::Timer T;
  T.start();
  source->setCheckQuality(checkQuality);
  source->transfer(target, method);
  T.stop();
  std::cout << "Time spent transferring data (ms) " << T.elapsed() << std::endl;
  target->write(ofname);
}


TransferDriver::TransferDriver(const std::string &srcmsh,
                               const std::string &trgmsh,
                               const std::string &method,
                               const std::vector<std::string> &arrayNames,
                               const std::string &ofname,
                               bool checkQuality)
{
  source = meshBase::Create(srcmsh);
  target = meshBase::Create(trgmsh);
  nemAux::Timer T;
  T.start();
  source->setCheckQuality(checkQuality);
  source->transfer(target, method, arrayNames);
  //source->write("new.vtu");
  T.stop();
  std::cout << "Time spent transferring data (ms) " << T.elapsed() << std::endl;
  target->write(ofname);
  std::cout << "TransferDriver created" << std::endl;
}


TransferDriver::~TransferDriver()
{
  delete source;
  delete target;
  std::cout << "TransferDriver destroyed" << std::endl;
}


TransferDriver *TransferDriver::readJSON(const jsoncons::json &inputjson)
{
  std::string srcmsh;
  std::string trgmsh;
  std::string outmsh;
  std::string method;
  std::string transferAll;
  std::string checkQual;
  bool transferall = true;
  bool checkQuality = false;
  std::vector<std::string> arrayNames;

  srcmsh = inputjson["Mesh File Options"]["Input Mesh Files"]["Source Mesh"].as<std::string>();
  trgmsh = inputjson["Mesh File Options"]["Input Mesh Files"]["Target Mesh"].as<std::string>();
  outmsh = inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  method = inputjson["Transfer Options"]["Method"].as<std::string>();

  transferAll = inputjson["Transfer Options"]["Transfer All Arrays"].as<std::string>();

  if (transferAll == "False" || transferAll == "false")
  {
    arrayNames = inputjson["Transfer Options"]["Array Names"].as<std::vector<std::string>>();
    transferall = false;
  }

  checkQual = inputjson["Transfer Options"]["Check Transfer Quality"].as<std::string>();
  if (checkQual == "True" || checkQual == "true")
  {
    checkQuality = true;
  }

  TransferDriver *trnsdrvobj;
  if (transferall)
  {
    trnsdrvobj = new TransferDriver(srcmsh, trgmsh, method,
                                    outmsh, checkQuality);
  }
  else
  {
    std::cout << "Transferring selected arrays:" << std::endl;
    for (const auto &arrayName : arrayNames)
    {
      std::cout << "\t" << arrayName << "\n";
    }
    trnsdrvobj = new TransferDriver(srcmsh, trgmsh, method,
                                    arrayNames,
                                    outmsh, checkQuality);
  }

  return trnsdrvobj;
}


TransferDriver *TransferDriver::readJSON(const std::string &ifname)
{
  std::ifstream inputStream(ifname);
  if (!inputStream.good() || nemAux::find_ext(ifname) != ".json")
  {
    std::cerr << "Error opening file " << ifname << std::endl;
    exit(1);
  }
  if (nemAux::find_ext(ifname) != ".json")
  {
    std::cerr << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  // checking if array
  if (inputjson.is_array())
  {
    std::cerr
        << "Warning: Input is an array. Only first element will be processed"
        << std::endl;
    return TransferDriver::readJSON(inputjson[0]);
  }
  else
  {
    return TransferDriver::readJSON(inputjson);
  }
}
