#ifndef MESHQUALITYDRIVER_H
#define MESHQUALITYDRIVER_H

#include "nemosys_export.h"
#include "NemDriver.H"

#include <string>

class NEMOSYS_EXPORT MeshQualityDriver : public NemDriver {
  public:
    MeshQualityDriver() : mesh(nullptr) {};
    MeshQualityDriver(const std::string &_mesh,
                      const std::string &ofname);
    ~MeshQualityDriver() override;

    static MeshQualityDriver *readJSON(const jsoncons::json &inputjson);

  private:
    meshBase *mesh;
};

#endif
