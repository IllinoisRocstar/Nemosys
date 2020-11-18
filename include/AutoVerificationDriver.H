#ifndef NEMOSYS_AUTOVERIFICATIONDRIVER_H
#define NEMOSYS_AUTOVERIFICATIONDRIVER_H

#include <string>

#include "nemosys_export.h"

#include "NemDriver.H"
#include "OrderOfAccuracy.H"
#include "meshBase.H"

#ifdef HAVE_OPENMP
#  include <omp.h>
#endif

namespace NEM {
namespace DRV {

  class NEMOSYS_EXPORT AutoVerificationDriver : public NemDriver {
  public:
    AutoVerificationDriver(meshBase *coarseMesh, meshBase *fineMesh,
                           meshBase *finerMesh, std::vector<int> arrayIDs,
                           std::string transferType, int numThreads = 2);

  public:
    static AutoVerificationDriver *readJSON(const std::string &ifname);

    static AutoVerificationDriver *readJSON(const jsoncons::json &inputjson);

  private:
    std::shared_ptr<OrderOfAccuracy> oac;
  };

#endif

} // namespace DRV
} // namespace NEM

