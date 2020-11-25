#ifndef NEMOSYS_OMEGAHREFINEDRIVER_H_
#define NEMOSYS_OMEGAHREFINEDRIVER_H_

#include "RefineDriver.H"
#include "geoMeshBase.H"

namespace NEM {
namespace DRV {

/**
 * @class omegahRefineDriver
 * @brief Driver to perform adaptive mesh refinement using Omega_h library
 */
class NEMOSYS_EXPORT omegahRefineDriver : public RefineDriver {
 public:
  omegahRefineDriver(NEM::MSH::geoMeshBase *in_mesh,
                     NEM::MSH::geoMeshBase *out_mesh,
                     const jsoncons::json &refine_opts);

  static omegahRefineDriver *readJSON(const jsoncons::json &inputjson);
};

}  // namespace DRV
}  // namespace NEM

#endif  // NEMOSYS_OMEGAHREFINEDRIVER_H_
