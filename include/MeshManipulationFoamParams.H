#ifndef NEMOSYS_MESHMANIPULATIONFOAMPARAMS_H_
#define NEMOSYS_MESHMANIPULATIONFOAMPARAMS_H_

#include "nemosys_export.h"

#include <string>

/** @brief MeshManipulationFoamParams defines basic parameters needed for
 *         utilities in MeshManipulationFoam class. These parameters are
 *         assigned a value during user input parsing through JSON in
 *         PackMeshDriver.
 */
class NEMOSYS_EXPORT MeshManipulationFoamParams {
  // --- constructors and destructors
 public:
  /** @brief Standard MeshManipulationFoamParams constructor
   */
  MeshManipulationFoamParams();

  /** @brief Standard MeshManipulationFoamParams destructor
   */
  ~MeshManipulationFoamParams() = default;

  // --- parameter variables
 public:
  // --- Booleans
  /** @brief Enables/Disables surfLambdaMuSmooth utility
             (Default is off)
   */
  bool _doSurfaceLMSmth;  // Enables SurfaceLambdaMuSmooth

  /** @brief Enables/Disables splitMshRegions utility
   *         (Default is off)
   */
  bool _doSplitMshRegs;  // Enables SplitMeshRegions

  /** @brief Enables/Disables mergeMeshes utility
   *         (Default is off)
   */
  bool _doMergeMsh;  // Enables MergeMeshes

  /** @brief Enables/Disables createPatch utility
   *         (Default is off)
   */
  bool _doCreatePtchs;  // Enables CreatePatch

  /** @brief Enables/Disables foamToSurface utility
   *         (Default is off)
   */
  bool _doFoam2Surf;  // Enables FoamToSurface

  /** @brief Enables/Disables surfSpltByTopology utility
   *         (Default is off)
   */
  bool _doSurfSplit;  // Enables SurfaceSplitByManifold

  // --- SurfaceLambdaMuSmooth
  /** @brief If enabled, allows adding feature file for surfLambdaMuSmooth
   *         (Default is off)
   */
  bool _addFeatureFile;  // Add feature file for smoothing

  /** @brief Input surface file for surfLambdaMuSmooth utility
   */
  std::string slmssurfaceFile;  // Defines input surface file

  /** @brief Output surface file for surfLambdaMuSmooth utility
   */
  std::string slmsoutputFile;  // Defines output surface file

  /** @brief Lambda (Range 0....1) for surfLambdaMuSmooth smoothing
   */
  double lambda;  // Lambda (0...1) for smoothing

  /** @brief Mu (Range 0....1) for surfLambdaMuSmooth smoothing.
   */
  double mu;  // Mu (0...1) for smoothing

  /** @brief Number of smoothing iterations for surfLambdaMuSmooth
   */
  int slmsIterations;  // Number of smoothing iterations

  // --- splitMeshRegions
  /** @brief Enables overwriting mesh in splitMshRegions
   *        (Default is on)
   */
  bool _overwriteMsh;  // Enable overwriting mesh

  /** @brief Enables using of cellZones option in splitMshRegions
   *        (Default is on)
   */
  bool _cellZones;  // Enable cellZones feature

  // --- mergeMeshes
  /** @brief Defines master mesh region for mergeMeshes
   */
  std::string masterCase;  // Defines master mesh region

  /** @brief Defines total number of domains (meshes) to merge using mergeMeshes
   */
  int numDomains;  // Defines number of disconnected domains

  /** @brief Defines name of first mesh region to add to master mesh
   */
  std::string addCase;  // Defines name of slave region to add

  /** @brief Enables overwriting merged mesh in mergeMeshes
   *        (Default is on)
   */
  bool _overwriteMergeMsh;  // Enable overwriting mesh

  /** @brief Specifies directory for master mesh
   */
  std::string masterCasePath;  // Defines directory of master region

  /** @brief Specifies directory for mesh to be merge with master mesh
   */
  std::string addCasePath;  // Defines directory of slave region

  // --- createPatch
  /** @brief Enables overwrite option for createPatch utility
   *        (Default is on)
   */
  bool _overwritecpMsh;  // Enable overwriting mesh

  /** @brief Defines surrounding patch name for createPatch utility
   *        (this is specific to "Pack Mesh" workflow)
   */
  std::string surroundingName;  // Patch name for surrounding

  /** @brief Defines packs patch name for createPatch utility
   *        (this is specific to "Pack Mesh" workflow)
   */
  std::string packsName;  // Patch name for packs

  /** @brief Defines surrounding patch type for createPatch utility
   *        (this is specific to "Pack Mesh" workflow)
   */
  std::string srrndngPatchType;  // Patch type - surrounding (wall, symm., etc.)

  /** @brief Defines packs patch type for createPatch utility
   *        (this is specific to "Pack Mesh" workflow)
   */
  std::string packsPatchType;  // Patch type - packs (wall, symm., etc.)

  /** @brief  Defines surrounding mesh path from main directory
   */
  std::string pathSurrounding;  // Path for surrounding mesh

  /** @brief Defines pack mesh path from main directory
   */
  std::string pathPacks;  // Path for pack mesh

  // --- foamToSurface
  /** @brief Output surface file name.extension for foamToSurface utility
   */
  std::string outSurfName;  // Output surface file name (Accepts path..
                            // ..from current directory)

  // --- surfaceSplitByTopology
  /** @brief Input surface file for surfSpltByTopology utility
   */
  std::string surfFile;  // Input surface file

  /** @brief Output surface file for surfSpltByTopology utility
   */
  std::string outSurfFile;  // Output surface file
};

#endif  // NEMOSYS_MESHMANIPULATIONFOAMPARAMS_H_
