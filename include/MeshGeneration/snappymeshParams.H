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
#ifndef NEMOSYS_SNAPPYMESHPARAMS_H_
#define NEMOSYS_SNAPPYMESHPARAMS_H_

#ifdef HAVE_CFMSH

#  include "nemosys_export.h"
#  include "MeshGeneration/meshingParams.H"

#  include <map>
#  include <memory>
#  include <string>
#  include <tuple>
#  include <vector>

/**
    @brief A structure to respresent regions
          refining capability in snappymeshGen.
**/
struct NEMOSYS_EXPORT shmRegionRefine  // Refinement Regions
{
  std::string patchName; /**< Patch name */
  std::string mode;      /**< Inside, Outside */
  double minLvl;         /**< Minimum refinement level */
  int maxLvl;            /**< Maximum refinement level */
};

/**
    @brief A structure to respresent shape patches
          defining capability in snappymeshGen.
**/
struct NEMOSYS_EXPORT shmSearchableShape {
  std::string patchName;
  virtual std::string getType() const = 0;
  virtual ~shmSearchableShape() = default;

 protected:
  shmSearchableShape() = default;
  shmSearchableShape(const shmSearchableShape &) = default;
  shmSearchableShape(shmSearchableShape &&) = default;
  shmSearchableShape &operator=(const shmSearchableShape &) = default;
  shmSearchableShape &operator=(shmSearchableShape &&) = default;
};

struct NEMOSYS_EXPORT shmSearchableBox : public shmSearchableShape {
  std::string getType() const override { return "searchableBox"; }
  std::array<double, 3> minBound;
  std::array<double, 3> maxBound;
};

struct NEMOSYS_EXPORT shmSearchableCylinder : public shmSearchableShape {
  std::string getType() const override { return "searchableCylinder"; }
  std::array<double, 3> axisPoint1;
  std::array<double, 3> axisPoint2;
  double radius;
};

struct NEMOSYS_EXPORT shmSearchableSphere : public shmSearchableShape {
  std::string getType() const override { return "searchableSphere"; }
  std::array<double, 3> center;
  double radius;
};

/**
    @brief A structure to represent geometry surface
          refining capability in snappymeshgGen
**/
struct NEMOSYS_EXPORT
    shmSurfRefine  // Surface Refinement (Add all shapes and their controls)
{
  std::string refPatchName;    /**< Refining region patch name */
  std::string patchType{"NO"}; /** < Patch Type in Mesh */
  double minLvl{1.};           /**< Minimum refinement level */
  int maxLvl;                  /**< Maximum refinement level */
};

/** @brief A structure for STL definition
 **/
struct NEMOSYS_EXPORT shmSTLDefinition {
  std::string snappyPatchName;
  std::string STLPatchName;
};

/** @brief A structure for feature edge refinement
 **/
struct NEMOSYS_EXPORT shmFeatureEdgeRef {
  std::string fileName;
  double minLvl{1.};
  int maxLvl;
};

/** @brief A struct for layer addition
 **/
struct NEMOSYS_EXPORT shmLayers {
  std::string patchName;
  int nSurfaceLayers{1};
  double expansionRatio{1.};
  double finalLayerThickness{1.};
  double firstLyrThickness{-1.};
  double thickness{-1.};
  double minThickness{1};
};

/** @brief Geometry definition
 **/
struct NEMOSYS_EXPORT shmGeo {
  /** @brief Boolean for if user wants whole STL as one solid of has
             different patches under STL
  **/
  bool withMultiPatches;

  /** @brief patch name for input surface file for when not @c _withMultiPatches
   **/
  std::string singleSolidPatch;

  /** @brief Vector for STL patches, if @c _withMultiPatches
   **/
  std::vector<shmSTLDefinition> stlPatchDefs;

  /** @brief Vector for searchable shapes
   **/
  std::vector<std::shared_ptr<shmSearchableShape>> srchShape;
};

/** @brief Castellated Mesh Controls
 **/
struct NEMOSYS_EXPORT shmCastMeshControls {
  /** @brief Enables defining cellzones and facezones
   **/
  bool withCellZones{false};  // Defining cellzones and facezones enabled?

  /** @brief Gap level increment
   **/
  int castMeshGpLvl{1};

  /** @brief max global cells allowed in mesh
   **/
  int maxLCells{2000000};  // max global cells

  /** @brief max local cells allowed (on 1 processor)
   **/
  int maxGCells{4000000};  // max local cells (on 1 processor)

  /** @brief minimum refinement cells
   **/
  int minRefCells{0};  // minimum refinement cells

  /** @brief number of cells between levels
   **/
  int cellsBetnLvls{3};  // number of cells between levels

  /** @brief minimum surface refinement
   **/
  int refSurfLvlMin{0};  // minimum surface refinement

  /** @brief maximum surface refinement
   **/
  int refSurfLvlMax{0};  // maximum surface refinement

  /** @brief resolve feature angle
   **/
  double featAngle{60.};  // resolve feature angle

  /** @brief Gap level increment
   **/
  int gapLvlInc{1};

  /** @brief Planar angle castallated mesh
   **/
  int planarAngle{30};

  /** @brief location in Mesh (cells are kept in regions accessinble from this
            location)
  **/
  std::array<double, 3> locMesh{0, 0, 0};

  /** @brief allows free standing zones if enabled
   **/
  bool alwFreeZone{true};  // allow free standing zones

  /** @brief Vector for feature edge refinement
   **/
  std::vector<shmFeatureEdgeRef> ftrEdge;

  /** @brief Vector for shmSurfRefine struct
   **/
  std::vector<shmSurfRefine> surfRefs;

  /** @brief Vector for shmRegionRefine struct
   **/
  std::vector<shmRegionRefine> geomRefs;
};

/** @brief Snapping Controls
 **/
struct NEMOSYS_EXPORT shmSnapControls {
  /** @brief Number of smoothing patches during snapping procedure
   **/
  int nSmoothPatch{4};  // nSmoothPatch

  /** @brief Snapping tolerance at shared interface (affects conformality of
            mesh)
  **/
  double tolerance{0.5};  // Snap Tolerance

  /** @brief Maximum iterations during snapping procedure before ending the
            process
  **/
  int nSolveIter{200};  // nSolveIter

  /** @brief Maximum number of relaxation iterations for snapping procedure
   **/
  int nRelaxIter{6};  // nRelaxIter

  /** @brief Feature snapping
   **/
  int nFeatureSnapIter{10};

  /** @brief Implicit feature snaping boolean
   **/
  bool implicitFeatureSnap{false};

  /** @brief Explicit feature snapping boolean
   **/
  bool explicitFeatureSnap{false};

  /** @brief Multi-region feature snapping boolean
   **/
  bool multiRegionFeatureSnap{false};
};

/** @brief Layer Controls
 **/
struct NEMOSYS_EXPORT shmLayerControls {
  /** @brief Enables relative sizes option during layer addition
   **/
  bool relSize{true};  // Relative Sizes

  /** @brief Expansion Ratio for layer addition
   **/
  double expRatio{1.3};  // Expansion Ratio

  /** @brief Final Layer Thickness for layer addition
   **/
  double finLThick{1.};  // Final Layer Thickness

  /** @brief Minimum Thickness for layer addition
   **/
  double minThick{0.1};  // Minimum Thickness

  /** @brief First Layer Thickness
   **/
  double firstLyrThickness{-1.};

  /** @brief Thickness
   **/
  double thickness{-1.};

  /** @brief Growth rate of successive layers
   **/
  int nGrow{0};

  /** @brief Feature Angle for layer addition
   **/
  double featAngle{30.};  // Feature Angle

  /** @brief number of relaxation interations for layer addition
   **/
  int relaxIter{3};  // number of relaxation interations

  /** @brief Number of smooth surface normals for layer addition
   **/
  int smthSurfNorm{1};  // # of smooth surface normals

  /** @brief Number of smooth normals for layer addition
   **/
  int smthNorm{3};  // # of smooth normals

  /** @brief Number of smooth thickness for layer addition
   **/
  int smthThick{2};  // # of smooth thickness

  /** @brief Maximum face thickness ratio for layer addition
   **/
  double maxFcTR{0.5};  // Maximum face thickness ratio

  /** @brief Maximum thickness to medial ratio for layer addition
   **/
  double maxThickTMR{1};  // Maximum thickness to medial ratio

  /** @brief Minimum medial axis angle for layer addition
   **/
  double minMedAngl{90.};  // Minimum medial axis angle

  /** @brief Number of buffer cells no extrude for layer addition
   **/
  int bufferCells{0};  // # of buffer cells no extrude

  /** @brief Number of layer interations for layer addition
   **/
  int nIter{50};  // # of layer interations

  /** @brief Relaxed iteration
   **/
  int nRelaxedIter{20};  // # of relaxed iteration

  /** @brief Slip feature angles
   **/
  int slipFeatureAngle{30};  // Slip feature angle

  /** @brief Limits # of steps walking away from the surface
   **/
  int nMedialAxisIter{-1};

  /** @brief Smooth displacement after medial axis determination
   **/
  int nSmoothDisplacement{-1};

  /** @brief Vector for Layers
   **/
  std::vector<shmLayers> layerVec;
};

/** @brief Mesh Quality Controls
 **/
struct NEMOSYS_EXPORT shmMeshQualityControls {
  /** @brief Maximum non-orthogonality
   **/
  int maxNonOrtho{65};  // Maximum non-orthogonality
  /** @brief Max Boundary Skewness
   **/
  double maxBndrySkew{20.};  // Max Boundary Skewness
  /** @brief Max Internal Skewness
   **/
  double maxIntSkew{4.};  // Max Internal Skewness
  /** @brief Max Concativity
   **/
  double maxConc{80.};  // Max Concativity
  /** @brief Minimum Cell Volume
   **/
  double minVol{1e-13};  // Minimum Cell Volume
  /** @brief Minimum Tet Quality
   **/
  double minTetQ{1e-15};  // Minimum Tet Quality
  /** @brief Minimum Area
   **/
  double minArea{-1.};  // Minimum Area
  /** @brief Minimum Twist
   **/
  double minTwist{0.02};  // Minimum Twist
  /** @brief Minimum Face Weight
   **/
  double minFaceW{0.05};  // Minimum Face Weight
  /** @brief Minimum Volume Ratio
   **/
  double minVolRto{0.01};  // Minimum Volume Ratio
  /** @brief Minimum Determinant
   **/
  double minDet{0.001};  // Minimum Determinant
  /** @brief Minimum Triangle Twist
   **/
  double minTriTwist{-1.};  // Minimum Triangle Twist
  /** @brief nSmoothScale
   **/
  int smoothScale{5};  // nSmoothScale
  /** @brief Error Reduction
   **/
  double errReduction{0.75};  // Error Reduction
};

/**
  @brief snappymeshParams contains all parameters essential for mesh generation
      using snappymeshGen class methods. These parameters are assigned a value
      at time of user input parsing through JSON in meshGemDriver.
**/
class NEMOSYS_EXPORT snappymeshParams : public meshingParams {
  // --- constructor and destructor
 public:
  /** @brief snappymeshParams standard constructor
   **/
  snappymeshParams() = default;

  /** @brief Geometry definition
   **/
  shmGeo geoDef;

  /** @brief Enables castellated mesh option for snappymeshGen
   **/
  bool withCastMesh{true};  // Boolean for castellated mesh

  /** @brief Castellated mesh controls
   **/
  shmCastMeshControls castMeshControls;

  /** @brief Enables surface snapping option for snappymeshGen
   **/
  bool withSnap{true};  // Boolean for snap

  /** @brief Snapping controls
   **/
  shmSnapControls snapControls;

  /** @brief Enables add layers option for snappymeshGen
   **/
  bool withLayers{false};  // Boolean for add layers

  /** @brief Layer controls
   **/
  shmLayerControls layerControls;

  /** @brief Input geometry STL name
   **/
  std::string geomFileName;  // Input geometry STL name

  /** @brief Mesh Quality Controls
   **/
  shmMeshQualityControls qualityControls;

  // --- Misc. General
  /** @brief merge tolerance for mesh
   **/
  double mergeTol{1e-6};  // Merge Tolerance

  /** @brief Boolean for packmesh
   **/
  bool isPackMesh{true};
};

#endif

#endif  // NEMOSYS_SNAPPYMESHPARAMS_H_
