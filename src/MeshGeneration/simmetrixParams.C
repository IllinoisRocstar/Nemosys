#include <simmetrixParams.H>

// constructor for simmetrixParams
simmetrixParams::simmetrixParams()
{
  meshSize = 0.05;
  anisoMeshCurv = 0.02;
  minCurvSize = 0.01;
  glbSizeGradRate = 0.1; 
  surfMshImprovGradRate = 0.05;
  surfMshImprovMinSize = 0.01;
} 
