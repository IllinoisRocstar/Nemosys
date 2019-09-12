#include "TransferDriver.H"

#include <iostream>
#include <string>

#include "AuxiliaryFunctions.H"

//----------------------- Transfer Driver ------------------------------------//
TransferDriver::TransferDriver(const std::string &srcmsh,
                               const std::string &trgmsh,
                               const std::string &method,
                               const std::string &ofname, bool checkQuality) {
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
                               const std::string &ofname, bool checkQuality) {
  source = meshBase::Create(srcmsh);
  target = meshBase::Create(trgmsh);
  std::cout << "TransferDriver created" << std::endl;

  nemAux::Timer T;
  T.start();
  source->setCheckQuality(checkQuality);
  source->transfer(target, method, arrayNames);
  // source->write("new.vtu");
  T.stop();

  std::cout << "Time spent transferring data (ms) " << T.elapsed() << std::endl;

  target->write(ofname);
}

TransferDriver::~TransferDriver() {
  delete source;
  delete target;
  std::cout << "TransferDriver destroyed" << std::endl;
}

TransferDriver *TransferDriver::readJSON(const jsoncons::json &inputjson) {
  std::string srcmsh =
      inputjson["Mesh File Options"]["Input Mesh Files"]["Source Mesh"]
          .as<std::string>();
  std::string trgmsh =
      inputjson["Mesh File Options"]["Input Mesh Files"]["Target Mesh"]
          .as<std::string>();
  std::string outmsh =
      inputjson["Mesh File Options"]["Output Mesh File"].as<std::string>();
  std::string method =
      inputjson["Transfer Options"]["Method"].as<std::string>();

  bool transferAll =
      inputjson["Transfer Options"]["Transfer All Arrays"].as<bool>();

  std::vector<std::string> arrayNames;
  if (!transferAll)
    arrayNames = inputjson["Transfer Options"]["Array Names"]
                     .as<std::vector<std::string>>();

  bool checkQuality =
      inputjson["Transfer Options"]["Check Transfer Quality"].as<bool>();

  TransferDriver *trnsdrvobj;
  if (transferAll) {
    trnsdrvobj =
        new TransferDriver(srcmsh, trgmsh, method, outmsh, checkQuality);
  } else {
    std::cout << "Transferring selected arrays:" << std::endl;
    for (const auto &arrayName : arrayNames)
      std::cout << "\t" << arrayName << "\n";
    trnsdrvobj = new TransferDriver(srcmsh, trgmsh, method, arrayNames, outmsh,
                                    checkQuality);
  }

  return trnsdrvobj;
}

TransferDriver *TransferDriver::readJSON(const std::string &ifname) {
  if (nemAux::find_ext(ifname) != ".json") {
    std::cerr << "Input File must be in .json format" << std::endl;
    exit(1);
  }

  std::ifstream inputStream(ifname);
  if (!inputStream.good()) {
    std::cerr << "Error opening file " << ifname << std::endl;
    exit(1);
  }

  jsoncons::json inputjson;
  inputStream >> inputjson;

  // checking if array
  if (inputjson.is_array()) {
    std::cerr
        << "Warning: Input is an array. Only first element will be processed\n";
    return TransferDriver::readJSON(inputjson[0]);
  } else {
    return TransferDriver::readJSON(inputjson);
  }
}
