#include "Drivers/InputGenDriver.H"

#include <iostream>
#include "AuxiliaryFunctions.H"
#include "ep16Post.H"
#include "ep16Prep.H"

namespace NEM {
namespace DRV {

InputGenDriver::InputGenDriver(std::string service, jsoncons::json opts)
    : service_(std::move(service)), opts_(std::move(opts)) {}

const std::string &InputGenDriver::getService() const { return service_; }

void InputGenDriver::setService(std::string service) {
  this->service_ = std::move(service);
}

const jsoncons::json &InputGenDriver::getOpts() const { return opts_; }

void InputGenDriver::setOpts(jsoncons::json opts) {
  this->opts_ = std::move(opts);
}

jsoncons::string_view InputGenDriver::getProgramType() const {
  return programType;
}

void InputGenDriver::execute() const {
  std::string srvName = this->service_;
  nemAux::toLower(srvName);
  nemAux::Timer T;
  T.start();
  if (srvName == "epic_2016") {
#ifdef HAVE_EPIC
    ep16Prep::readJSON(this->opts_);
#else
    std::cerr << "Compile the code with EPIC module enabled." << std::endl;
    exit(0);
#endif
  } else if (srvName == "epic_2016_post") {
#ifdef HAVE_EPIC
    int ret;
    NEM::EPC::ep16Post::readJSON(this->opts_, ret);
#else
    std::cerr << "Compile the code with EPIC module enabled." << std::endl;
    exit(0);
#endif
  } else {
    std::cerr << "The input generation service " << srvName << "is unsupported."
              << std::endl;
  }
  T.stop();
}

}  // namespace DRV
}  // namespace NEM
