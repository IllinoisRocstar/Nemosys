#ifndef NEMOSYS_CFMESHQUALITYPARAMS_H_
#define NEMOSYS_CFMESHQUALITYPARAMS_H_

#ifdef HAVE_CFMSH

#include "nemosys_export.h"

#include <string>

#include "meshQualityParams.H"

class NEMOSYS_EXPORT cfmshQualityParams : public meshQualityParams {
 public:
  // initializing params with defaults
  cfmshQualityParams();
  ~cfmshQualityParams() override = default;

  // meshOptimizer
 public:
  int nIterations;          // number of optimization iterations
  int nLoops;               // number of inner loops in optimization
  double qualThrsh;         // quality threshold
  int nSrfItr;              // number of surface iterations
  std::string consCellSet;  // name of constrained cellSet
  bool _withConstraint;
};

#endif  // HAVE_CFMESH

#endif  // NEMOSYS_CFMESHQUALITYPARAMS_H_
