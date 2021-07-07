#ifndef NEMOSYS_CFMESHQUALITYPARAMS_H_
#define NEMOSYS_CFMESHQUALITYPARAMS_H_

#ifdef HAVE_CFMSH

#include "nemosys_export.h"
#include "MeshQuality/meshQualityParams.H"

#include <jsoncons/config/jsoncons_config.hpp>  // for jsoncons::optional
#include <jsoncons/json_traits_macros.hpp>
#include <string>

class NEMOSYS_EXPORT cfmeshQualityParams : public meshQualityParams {
 public:
  // initializing params with defaults
  cfmeshQualityParams() = default;

  int nIterations{50};          // number of optimization iterations
  int nLoops{10};               // number of inner loops in optimization
  double qualThrsh{0.1};         // quality threshold
  int nSrfItr{2};              // number of surface iterations
  jsoncons::optional<std::string> consCellSet{};  // name of constrained cellSet

  JSONCONS_TYPE_TRAITS_FRIEND

 protected:
  static constexpr const char *method = "meshOptimizer";
};

#endif  // HAVE_CFMSH

#endif  // NEMOSYS_CFMESHQUALITYPARAMS_H_
