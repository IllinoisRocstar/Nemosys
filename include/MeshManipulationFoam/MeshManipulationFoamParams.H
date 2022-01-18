/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_MESHMANIPULATIONFOAMPARAMS_H_
#define NEMOSYS_MESHMANIPULATIONFOAMPARAMS_H_

#include "nemosys_export.h"

#include <string>
#include <vector>

/** @brief MeshManipulationFoamParams defines basic parameters needed for
 *         utilities in MeshManipulationFoam class. These parameters are
 *         assigned a value during user input parsing through JSON in
 *         PackMeshDriver.
 */
class NEMOSYS_EXPORT MeshManipulationFoamParams {
 public:
  struct NEMOSYS_EXPORT SurfaceLambdaMuSmooth {
    // --- SurfaceLambdaMuSmooth
    /** @brief If enabled, allows adding feature file for surfLambdaMuSmooth
     *         (Default is off)
     */
    bool addFeatureFile;  // Add feature file for smoothing

    /** @brief Input surface file for surfLambdaMuSmooth utility
     */
    std::string slmssurfaceFile;  // Defines input surface file

    /** @brief Output surface file for surfLambdaMuSmooth utility
     */
    std::string slmsoutputFile;  // Defines output surface file

    /** @brief Lambda (Range 0....1) for surfLambdaMuSmooth smoothing
     */
    double lambda_{1.};  // Lambda (0...1) for smoothing

    /** @brief Mu (Range 0....1) for surfLambdaMuSmooth smoothing.
     */
    double mu{1.};  // Mu (0...1) for smoothing

    /** @brief Number of smoothing iterations for surfLambdaMuSmooth
     */
    int slmsIterations{50};  // Number of smoothing iterations
  };

  struct NEMOSYS_EXPORT SplitMeshRegions {
    // --- splitMeshRegions
    /** @brief Enables overwriting mesh in splitMshRegions
     *        (Default is on)
     */
    bool overwriteMsh{true};  // Enable overwriting mesh

    /** @brief Enables using of cellZones option in splitMshRegions
     *        (Default is on)
     */
    bool cellZones{true};  // Enable cellZones feature
  };

  struct NEMOSYS_EXPORT MergeMeshes {
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
    bool overwriteMergeMsh{true};  // Enable overwriting mesh

    /** @brief Specifies directory for master mesh
     */
    std::string masterCasePath;  // Defines directory of master region

    /** @brief Specifies directory for mesh to be merge with master mesh
     */
    std::string addCasePath;  // Defines directory of slave region
  };

  struct NEMOSYS_EXPORT CreatePatch {
    // --- createPatch
    /** @brief Enables overwrite option for createPatch utility
     *        (Default is on)
     */
    bool overwritecpMsh{true};  // Enable overwriting mesh

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
  };

  struct NEMOSYS_EXPORT FoamToSurface {
    // --- foamToSurface
    /** @brief Output surface file name.extension for foamToSurface utility
     */
    std::string outSurfName;  // Output surface file name (Accepts path from
                              // current directory)
  };

  struct NEMOSYS_EXPORT SurfaceSplitByManifold {
    // --- surfaceSplitByTopology
    /** @brief Input surface file for surfSpltByTopology utility
     */
    std::string surfFile;  // Input surface file

    /** @brief Output surface file for surfSpltByTopology utility
     */
    std::string outSurfFile;  // Output surface file

    /** @brief A vector containing all region names for packs
     **/
    std::vector<std::string> pckRegionNames;
  };

  SurfaceLambdaMuSmooth surfLMSmoothParams;

  SplitMeshRegions splitMeshRegParams;

  MergeMeshes mergeMeshesParams;

  CreatePatch createPatchParams;

  FoamToSurface foamToSurfParams;

  SurfaceSplitByManifold surfSplitParams;
};

#endif  // NEMOSYS_MESHMANIPULATIONFOAMPARAMS_H_
