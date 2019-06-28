#ifndef MESHGENDRIVER_H
#define MESHGENDRIVER_H

#include "nemosys_export.h"
#include "NemDriver.H"

#include "meshingParams.H"
#include "meshBase.H"

#include <string>
#include <memory>

class NEMOSYS_EXPORT MeshGenDriver : public NemDriver {
  public:
    MeshGenDriver(const std::string &ifname,
                  const std::string &meshEngine,
                  meshingParams *params,
                  const std::string &ofname);
    ~MeshGenDriver() override;

    static MeshGenDriver *readJSON(const jsoncons::json &inputjson);
    static MeshGenDriver *readJSON(const std::string &ifname,
                                   const std::string &ofname,
                                   const jsoncons::json &inputjson);
    std::shared_ptr<meshBase> getNewMesh() const;

  private:
    std::shared_ptr<meshBase> mesh;
    meshingParams *params;
};

#endif
