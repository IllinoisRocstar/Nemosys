#ifndef CONVERSIONDRIVER_H
#define CONVERSIONDRIVER_H

#include "nemosys_export.h"
#include "NemDriver.H"

#ifdef HAVE_EXODUSII
  #include "exoMesh.H"
#endif

#include <string>

class NEMOSYS_EXPORT ConversionDriver : public NemDriver {
  public:
    ConversionDriver(const std::string &srcmsh,
                     const std::string &trgmsh,
                     const std::string &method,
                     const std::string &ofname,
                     const jsoncons::json &inputjson);
    ~ConversionDriver() override;

    static ConversionDriver *readJSON(const jsoncons::json &inputjson);
    static ConversionDriver *readJSON(const std::string &ifname);

#ifdef HAVE_EXODUSII
  private:
    void genExo(const jsoncons::json &opts,
                const std::string &fname);
    void procExo(const jsoncons::json &ppJson,
                 const std::string &fname,
                 EXOMesh::exoMesh *em);
#endif

  private:
    meshBase *source;
};

#endif
