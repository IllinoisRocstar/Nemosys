#ifndef MESHMANIPULATIONFOAMPARAMS_H
#define MESHMANIPULATIONFOAMPARAMS_H

#include <string>
#include <vector>

class MeshManipulationFoamParams
{
  public:
  MeshManipulationFoamParams();
  ~MeshManipulationFoamParams(){};

  public:
  //Booleans
  bool _doSurfaceLMSmth;      // Enables SurfaceLamndaMuSmooth
  bool _doSplitMshRegs;     // Enables SplitMeshRegions
  bool _doMergeMsh;       // Enables MergeMeshes
  bool _doCreatePtchs;      // Enables CreatePatch
  bool _doFoam2Surf;        // Enables FoamToSurface
  bool _doSurfSplit;        // Enables SurfaceSplitByManifold

  // SurfaceLambdaMuSmooth
  bool _addFeatureFile;     // Add feature file for smoothing
  std::string slmssurfaceFile;  // Defines input surface file
  std::string slmsoutputFile;   // Defines output surface file
  double lambda;          // Lambda (0...1) for smoothing
  double mu;            // Mu (0...1) for smoothing
  int slmsIterations;       // Number of smoothing iterations

  //splitMeshRegions
  bool _overwriteMsh;       // Enable overwriting mesh
  bool _cellZones;        // Enable cellZones feature

  //mergeMeshes
  std::string masterCase;     // Defines master mesh region
  int numDomains;         // Defines number of disconnected domains
  std::string addCase;      // Defines name of slave region to add
  bool _overwriteMergeMsh;    // Enable overwriting mesh
  std::string masterCasePath;   // Defines directory of master region
  std::string addCasePath;    // Defines directory of slave region

  //createPatch
  bool _overwritecpMsh;     // Enable overwriting mesh
  std::string surroundingName;  // Patch name for surrounding
  std::string packsName;      // Patch name for packs
  std::string srrndngPatchType; // Patch type - surrounding(wall,symm.,etc)
  std::string packsPatchType;   // Patch type - packs(wall,symm.,etc)
  std::string pathSurrounding;  // Path for surrounding mesh
  std::string pathPacks;      // Path for pack mesh

  //foamToSurface
  std::string outSurfName;    // Output surface file name (Accepts path..
                  // ..from current directory)

  //surfaceSplitByTopology
  std::string surfFile;     // Input surface file
  std::string outSurfFile;    // Output surface file
};

#endif