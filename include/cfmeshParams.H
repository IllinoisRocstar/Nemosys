#ifndef CFMESHPARAMS_H
#define CFMESHPARAMS_H

#ifdef HAVE_CFMSH

#include "meshingParams.H"

#include <string>
#include <vector>
#include <map>
#include <tuple>

struct cfmObjRef
{
  std::string name;
  std::map<std::string, std::string> params; // all other parameters
};

struct cfmLclRefPatch
{
  std::string patchName;    // patch name
  double cellSize;          // patch cell size
  int aditRefLvls;          // additional refinement levels
  double refThickness;      // local refinement thickness
};

struct cfmPtchBndLyr
{
  std::string patchName;    // patch name
  int blNLyr;               // number of boundary layers
  double blThkRto;          // thickness ratio
  double maxFrstLyrThk;     // thickness of first boundary layer
  bool alwDiscont;          // allowing for discontinuity can be used for patches
};

typedef std::tuple<std::string, std::string, std::string> cfmNewPatch;
struct cfmRenBndry
{
  std::string defName;
  std::string defType;
  std::vector<cfmNewPatch> newPatches;
};

class cfmeshParams : public meshingParams
{
  public:
    // initializing params with defaults
    cfmeshParams();

    ~cfmeshParams() override = default;

    // general
    bool _alwDiscDomains;       // Multi-region meshing
    std::string geomFilePath;   // name of geometry file and full path
    std::string generator;      // mesh generation method to be used
    double maxCellSize;         // maximum cell size
    double minCellSize;         // minimum cell size
    double bndryCellSize;       // boundary cell size
    double bndryCellSizeRefThk; // boundary cell size refinement thickness
    bool keepCellIB;            // keep cells intersecting boundary
    bool chkGluMsh;             // check for glued mesh

    // boundary layer
    bool _withBndLyr;
    bool _withBndLyrPtch;
    int blNLyr;                 // number of boundary layers
    double blThkRto;            // thickness ratio
    double maxFrstLyrThk;       // thickness of first boundary layer
    bool alwDiscont;            // allowing for discontinuity can be used for patches

    // surfaceFeatureEdges
    bool _withSrfEdg;
    double srfEdgAng;           // angle used for feature extraction

    // object refinements
    bool _withObjRfn;
    std::string refTyp;         // type of refinement size field
    double refCellSize;         // target cell size
    std::string refP0;          // P0 point coordinates eg (0,0,0)
    std::string refP1;          // P1 point coordinates eg (1,1,1)
    double refRad0;             // radius 0
    double refRad1;             // radius 1

    // mesh quality improvement
    bool _withMshQlt;
    int qltNItr;                // number of iterations
    int qltNLop;                // number of loops
    double qltQltThr;           // quality threshold
    int qltNSrfItr;             // number of surface iterations
    std::string qltConCelSet;   // name of constrained cellSet

    // local refinement 
    bool _withLclRef;
    std::vector<cfmObjRef> objRefLst;         // objects to refine
    std::vector<cfmPtchBndLyr> blPatches;     // define additional constrains
                                              // on boundary layers
    std::vector<cfmLclRefPatch> refPatches;   // define additional constrains
                                              // on refinement of given patches
    // rename boundaries
    bool _withRenBndry;
    cfmRenBndry renBndry;                     // rename boundaries
};

#endif

#endif
