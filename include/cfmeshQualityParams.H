#ifndef CFMESHQUALITYPARAMS_H
#define CFMESHQUALITYPARAMS_H

#ifdef HAVE_CFMSH

#include "meshQualityParams.H"

#include <string>


class cfmshQualityParams : public meshQualityParams
{
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

#endif

#endif
