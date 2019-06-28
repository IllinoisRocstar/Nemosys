#include "cfmeshParams.H"


cfmeshParams::cfmeshParams()
{
  // optional values are defaulted to negative values
  maxCellSize = -1.;
  minCellSize = -1.;
  bndryCellSize = -1.;
  maxFrstLyrThk = -1.;
  alwDiscont = false;
  keepCellIB = -1;
  chkGluMsh = -1;
  blNLyr = 1;
  blThkRto = 1.0;
  srfEdgAng = 45.0;

  _withBndLyr = false;
  _withSrfEdg = false;
  _withObjRfn = false;
  _withObjRfn = false;
  _withMshQlt = false;
  _withLclRef = false;
  _withBndLyrPtch = false;
  _withRenBndry = false;
}
