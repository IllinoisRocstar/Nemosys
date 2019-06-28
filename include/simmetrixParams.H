#ifndef SIMMETRIXPARAMS_H
#define SIMMETRIXPARAMS_H

#include "nemosys_export.h"

#include "meshingParams.H"

#include <string>

/* possible feature options are :
    - geomsim_core
    - meshsim_surface
    - meshsim_volume
    - geomsim_discrete
    - meshsim_adv
    - meshsim_adapt */

class NEMOSYS_EXPORT simmetrixParams : public meshingParams
{
  public:
    simmetrixParams();

    ~simmetrixParams() override = default;

  public:
    std::string logFName;
    std::string features;
    std::string licFName;
    double meshSize;
    double anisoMeshCurv;
    double minCurvSize;
    double glbSizeGradRate;
    double surfMshImprovGradRate;
    double surfMshImprovMinSize;
};

#endif
