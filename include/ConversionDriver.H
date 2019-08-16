#ifndef NEMOSYS_CONVERSIONDRIVER_H_
#define NEMOSYS_CONVERSIONDRIVER_H_

#include "nemosys_export.h"

#include <string>

#include "NemDriver.H"
#ifdef HAVE_EXODUSII
#  include "exoMesh.H"
#endif

class NEMOSYS_EXPORT ConversionDriver : public NemDriver {
 public:
  ConversionDriver(const std::string &srcmsh, const std::string &trgmsh,
                   const std::string &method, const std::string &ofname,
                   const jsoncons::json &inputjson);
  ~ConversionDriver() override;

  static ConversionDriver *readJSON(const jsoncons::json &inputjson);
  static ConversionDriver *readJSON(const std::string &ifname);

#ifdef HAVE_EXODUSII
 private:
  static void genExo(const jsoncons::json &opts, const std::string &fname);
  static void procExo(const jsoncons::json &ppJson, const std::string &fname,
                      EXOMesh::exoMesh *em);
#endif

 private:
  meshBase *source;
};

#endif  // NEMOSYS_CONVERSIONDRIVER_H_
