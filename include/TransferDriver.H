#ifndef TRANSFERDRIVER_H
#define TRANSFERDRIVER_H

#include "nemosys_export.h"
#include "NemDriver.H"

class NEMOSYS_EXPORT TransferDriver : public NemDriver {
  public:
    TransferDriver(const std::string &srcmsh,
                   const std::string &trgmsh,
                   const std::string &method,
                   const std::string &ofname,
                   bool checkQuality);

    TransferDriver(const std::string &srcmsh,
                   const std::string &trgmsh,
                   const std::string &method,
                   const std::vector<std::string> &arrayNames,
                   const std::string &ofname,
                   bool checkQuality);

    static TransferDriver *readJSON(const jsoncons::json &inputjson);
    static TransferDriver *readJSON(const std::string &ifname);

    ~TransferDriver() override;

  private:
    meshBase *source;
    meshBase *target;
};

#endif
