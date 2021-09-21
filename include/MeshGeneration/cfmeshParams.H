#ifndef NEMOSYS_CFMESHPARAMS_H_
#define NEMOSYS_CFMESHPARAMS_H_

#ifdef HAVE_CFMSH

#include "nemosys_export.h"
#include "MeshGeneration/meshingParams.H"

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <jsoncons/config/jsoncons_config.hpp>  // for jsoncons::optional

struct NEMOSYS_EXPORT cfmObjRef {
  std::string name;
  std::map<std::string, std::string> params;  // all other parameters
};

struct NEMOSYS_EXPORT cfmLclRefPatch {
  std::string patchName;     // patch name
  double cellSize{-1.};      // patch cell size
  int aditRefLvls{-1};       // additional refinement levels
  double refThickness{-1.};  // local refinement thickness
};

struct NEMOSYS_EXPORT cfmPtchBndLyr {
  std::string patchName;      // patch name
  int blNLyr{-1};             // number of boundary layers
  double blThkRto{1.};        // thickness ratio
  double maxFrstLyrThk{-1.};  // thickness of first boundary layer
  bool alwDiscont{false};  // allowing for discontinuity can be used for patches
};

struct NEMOSYS_EXPORT cfmNewPatch {
  std::string name;
  std::string newName;
  std::string newType;
};

struct NEMOSYS_EXPORT cfmRenBndry {
  std::string defName;
  std::string defType;
  std::vector<cfmNewPatch> newPatches;
};

struct NEMOSYS_EXPORT cfmBoundaryLayer {
  int blNLyr{1};              // number of boundary layers
  double blThkRto{1.};        // thickness ratio
  double maxFrstLyrThk{-1.};  // thickness of first boundary layer
  // allowing for discontinuity can be used for patches
  bool alwDiscont{false};
  // define additional constrains on boundary layers
  std::vector<cfmPtchBndLyr> blPatches{};
};

struct NEMOSYS_EXPORT cfmMeshQual {
  int qltNItr;                       // number of iterations
  int qltNLop;                       // number of loops
  double qltQltThr;                  // quality threshold
  int qltNSrfItr;                    // number of surface iterations
  std::string qltConCelSet{"none"};  // name of constrained cellSet
};

struct NEMOSYS_EXPORT cfmSrfFeatEdge {
  double srfEdgAng{45.};
};

class NEMOSYS_EXPORT cfmeshParams : public meshingParams {
 public:
  // initializing params with defaults
  cfmeshParams() = default;

  // general
  bool alwDiscDomains{false};       // Multi-region meshing
  std::string geomFilePath;         // name of geometry file and full path
  std::string generator;            // mesh generation method to be used
  double maxCellSize{-1.};          // maximum cell size
  double minCellSize{-1.};          // minimum cell size
  double bndryCellSize{-1.};        // boundary cell size
  double bndryCellSizeRefThk{-1.};  // boundary cell size refinement thickness
  bool keepCellIB{false};           // keep cells intersecting boundary
  bool chkGluMsh{false};            // check for glued mesh

  // boundary layer
  jsoncons::optional<cfmBoundaryLayer> boundaryLayers;

  // surfaceFeatureEdges
  jsoncons::optional<cfmSrfFeatEdge> srfEdge{};

  // object refinements
  std::vector<cfmObjRef> objRefLst;  // objects to refine
  std::string refTyp;                // type of refinement size field
  double refCellSize;                // target cell size
  std::string refP0;                 // P0 point coordinates eg (0,0,0)
  std::string refP1;                 // P1 point coordinates eg (1,1,1)
  double refRad0;                    // radius 0
  double refRad1;                    // radius 1

  // mesh quality improvement
  jsoncons::optional<cfmMeshQual> improveMeshQuality{};

  // local refinement
  std::vector<cfmLclRefPatch> refPatches;  // define additional constrains
                                           // on refinement of given patches
  // rename boundaries
  jsoncons::optional<cfmRenBndry> renBndry;  // rename boundaries

  // Boolean for if the operation is packmesh
  bool isPackMesh{true};
};

#endif

#endif  // NEMOSYS_CFMESHPARAMS_H_
