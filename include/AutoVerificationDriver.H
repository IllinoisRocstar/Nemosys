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

/**
 * This class drives the automatic verification (AV) in Nemosys. Three
 * successively refined meshes of a given geometry are required to evaluate
 * whether the coarsest mesh is within the asymptotic range. The grid
 * convergence index (GCI) of each component of each selected field is
 * evaluated: if they are all approximately equal to 1 then the coarse mesh is
 * sufficiently refined for an accurate solution (and is within asymptotic
 * range).
 */
class NEMOSYS_EXPORT AutoVerificationDriver : public NemDriver {
 public:
  /**
   * @brief Creates AV driver given successively refined meshes of a
   * given geometry.
   * @param coarseMesh coarsest mesh
   * @param finerMesh fine mesh
   * @param finestMesh finest mesh
   * @param arrayIDs field ids for which we will evaluate the grid convergence
   * index (GCI)
   * @param transferType default is 'Consistent Interpolation',
   * alternatives include 'Conservative Surface Transfer' (if enabled, requires
   * IMPACT) and 'Conservative Volume Transfer' (if enabled)
   * @param targetGCI the target GCI value, set to 1.1 by default
   * @param numThreads sets number of threads used in transfer (OpenMP), set to
   * 2 by default
   */
  AutoVerificationDriver(meshBase *coarseMesh, meshBase *finerMesh,
                         meshBase *finestMesh, std::vector<int> arrayIDs,
                         std::string transferType, double targetGCI = 1.1,
                         int numThreads = 2);

 public:
  static AutoVerificationDriver *readJSON(const std::string &ifname);

  static AutoVerificationDriver *readJSON(const jsoncons::json &inputjson);

 private:
  std::shared_ptr<OrderOfAccuracy> oac;
};

}  // namespace DRV
}  // namespace NEM

#endif
