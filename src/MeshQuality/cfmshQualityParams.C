#include "cfmeshQualityParams.H"


cfmshQualityParams::cfmshQualityParams()
{
  // default values
  nIterations = 50;
  nLoops = 10;
  qualThrsh = 0.1;
  nSrfItr = 2;

  _withConstraint = false;
}
