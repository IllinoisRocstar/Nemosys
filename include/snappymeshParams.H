#ifndef SNAPPYMESHPARAMS_H
#define SNAPPYMESHPARAMS_H

#ifdef HAVE_CFMSH

#  include <map>
#  include <string>
#  include <tuple>
#  include <vector>
#  include "meshingParams.H"

/**
    @brief A structure to respresent regions
          refining capability in snappymeshGen.
**/
struct shmRegionRefine  // Refinement Regions
{
  std::string patchNm; /**< Patch name */
  std::string mode;    /**< Inside, Outside */
  double minLvl;       /**< Minimum refinement level */
  int maxLvl;          /**< Maximum refinement level */
};

/**
    @brief A structure to respresent shape patches
          defining capability in snappymeshGen.
**/
struct shmSearchableShapes {
  std::string patchNm;          /**< Patch name */
  std::string searchableName;   /**< Shape name */
  std::string shapeParameters1; /**< Points, Center etc */
  std::string shapeParameters2; /**< Points, Center etc */
  double rad;                   /**< Radius */
};

/**
    @brief A structure to represent geometry surface
          refining capability in snappymeshgGen
**/
struct shmSurfRefine  // Surface Refinement (Add all shapes and their controls)
{
  std::string refPatchNm; /**< Refining region patch name */
  std::string patchType;  /** < Patch Type in Mesh */
  double minLvl;          /**< Minimum refinement level */
  int maxLvl;             /**< Maximum refinement level */
};

/** @brief A structure for STL definition
 **/
struct shmSTLDefinition {
  std::string snappyPatchName;
  std::string STLPatchName;
};

/** @brief A structure for feature edge refinement
 **/
struct shmFeatureEdgeRef {
  std::string fileName;
  double minLvl;
  int maxLvl;
};

/** @brief A struct for layer addition
 **/
struct shmLayers {
  std::string patchName;
  int nSurfaceLayers;
  double expansionRatio;
  double finalLayerThickness;
  double firstLyrThickness;
  double thickness;
  double minThickness;
};

/**
  @brief snappymeshParams contains all parameters essential for mesh generation
      using snappymeshGen class methods. These parameters are assigned a value
      at time of user input parsing through JSON in meshGemDriver.
**/
class snappymeshParams : public meshingParams {
  // --- constructor and destructor
 public:
  /** @brief snappymeshParams standard constructor
   **/
  snappymeshParams();

  /** @brief snappymeshParams standard destructor
   **/
  ~snappymeshParams(){};

  // --- Booleans
  /** @brief Enables castellated mesh option for snappymeshGen
   **/
  bool _withCastMesh;  // Boolean for castellated mesh

  /** @brief Enables surface snapping option for snappymeshGen
   **/
  bool _withSnap;  // Boolean for snap

  /** @brief Enables add layers option for snappymeshGen
   **/
  bool _withLayers;  // Boolean for add layers

  /** @brief Enables user defined geometry refinement regions (searchable)
   **/
  bool _withGeomRefReg;  // Geometry refinement regions (searchable)

  /** @brief Enables defining cellzones and facezones
   **/
  bool _withCellZones;  // Defining cellzones and facezones enabled?

  /** @brief Enables user defined surface refinement regions
   **/
  bool _withSurfRefReg;  // Surface refinement regions

  /** @brief Boolean for if user wants whole STL as one solid of has
             different patches under STL
  **/
  bool _withMultiPatches;

  /** @brief Boolean for feature edge refinement
   **/
  bool _withFeatureEdgeFile;

  // --- Castellated Mesh Controls
  /** @brief Input geometry STL name
   **/
  std::string geomFileName;  // Input geometry STL name

  /** @brief patch name for input surface file
   **/
  std::string singleSolidPatch;  // patch name for same geometry file

  /** @brief optional non-default tolerance on intersections
   **/
  double geomTolerance;

  /** @brief optional depth of octree. Decrease only in case
   **/
  int maxTreeDepth;

  /** @brief surface refinement for whole geometry
   **/
  int wholeGeomMinRef;

  /** @brief surface refinement for whole geometry
   **/
  int wholeGeomMaxRef;

  /** @brief max global cells allowed in mesh
   **/
  int maxLCells;  // max global cells

  /** @brief max local cells allowed (on 1 processor)
   **/
  int maxGCells;  // max local cells (on 1 processor)

  /** @brief minimum refinement cells
   **/
  int minRefCells;  // minimum refinement cells

  /** @brief number of cells between levels
   **/
  int cellsBetnLvls;  // number of cells between levels

  /** @brief minimum surface refinement
   **/
  int refSurfLvlMin;  // minimum surface refinement

  /** @brief maximum surface refinement
   **/
  int refSurfLvlMax;  // maximum surface refinement

  /** @brief resolve feature angle
   **/
  double featAngle;  // resolve feature angle

  /** @brief location in Mesh (cells are kept in regions accessinble from this
            location)
  **/
  double locMeshX;  // location in Mesh

  /** @brief location in Mesh (cells are kept in regions accessinble from this
            location)
  **/
  double locMeshY;  // location in Mesh

  /** @brief location in Mesh (cells are kept in regions accessinble from this
            location)
  **/
  double locMeshZ;  // location in Mesh

  /** @brief allows free standing zones if enabled
   **/
  bool _alwFreeZone;  // allow free standing zones

  /** @brief Gap level increment
   **/
  int gPLvlInc;

  /** @brief Gap level increment
   **/
  int castMeshGpLvl;

  /** @brief Planar angle castallated mesh
   **/
  int planarAngle;

  // --- Snap Controls
  /** @brief Number of smoothing patches during snapping procedure
   **/
  int snapSmthPatch;  // nSmoothPatch

  /** @brief Snapping tolerance at shared interface (affects conformality of
            mesh)
  **/
  double snapTol;  // Snap Tolerance

  /** @brief Maximum iterations during snapping procedure before ending the
            process
  **/
  int solveSnapIter;  // nSolveIter

  /** @brief Maximum number of relaxation iterations for snapping procedure
   **/
  int relaxSnapIter;  // nRelaxIter

  /** @brief Feature snapping
   **/
  int nFeatureSnapIter;

  /** @brief Implicit feature snaping boolean
   **/
  bool implicitFeatureSnap;

  /** @brief Explicit feature snapping boolean
   **/
  bool explicitFeatureSnap;

  /** @brief Multi-region feature snapping boolean
   **/
  bool multiRegionFeatureSnap;

  // --- Layer Controls
  /** @brief Enables relative sizes option during layer addition
   **/
  bool _relSize;  // Relative Sizes

  /** @brief Expansion Ratio for layer addition
   **/
  double expRatio;  // Expansion Ratio

  /** @brief Final Layer Thickness for layer addition
   **/
  double finLThick;  // Final Layer Thickness

  /** @brief Minimum Thickness for layer addition
   **/
  double minThick;  // Minimum Thickness

  /** @brief First Layer Thickness
   **/
  double firstLyrThickness;

  /** @brief Thickness
   **/
  double thickness;

  /** @brief Growth rate of successive layers
   **/
  int nGrow;

  /** @brief Feature Angle for layer addition
   **/
  double lyrFeatAngle;  // Feature Angle

  /** @brief number of relaxation interations for layer addition
   **/
  int lyrRelaxIter;  // number of relaxation interations

  /** @brief Number of smooth surface normals for layer addition
   **/
  int lyrSmthSurfNorm;  // # of smooth surface normals

  /** @brief Number of smooth normals for layer addition
   **/
  int lyrSmthNorm;  // # of smooth normals

  /** @brief Number of smooth thickness for layer addition
   **/
  int lyrSmthThick;  // # of smooth thickness

  /** @brief Maximum face thickness ratio for layer addition
   **/
  double lyrMaxFcTR;  // Maximum face thickness ratio

  /** @brief Maximum thickness to medial ratio for layer addition
   **/
  double lyrMaxThickTMR;  // Maximum thickness to medial ratio

  /** @brief Minimum medial axis angle for layer addition
   **/
  double lyrMinMedAngl;  // Minimum medial axis angle

  /** @brief Number of buffer cells no extrude for layer addition
   **/
  int lyrBuffrCells;  // # of buffer cells no extrude

  /** @brief Number of layer interations for layer addition
   **/
  int lyrIter;  // # of layer interations

  /** @brief Relaxed iteration
   **/
  int nRelaxedIter;  // # of relaxed iteration

  /** @brief Slip feature angles
   **/
  int slipFeatureAngle;  // Slip feature angle

  /** @brief Limits # of steps walking away from the surface
   **/
  int nMedialAxisIter;

  /** @brief Smooth displacement after medial axis determination
   **/
  int nSmoothDisplacement;

  // --- Mesh Quality Controls
  /** @brief Maximum non-orthogonality
   **/
  int qcMaxNOrtho;  // Maximum non-orthogonality
  /** @brief Max Boundary Skewness
   **/
  double qcMaxBndrySkew;  // Max Boundary Skewness
  /** @brief Max Internal Skewness
   **/
  double qcMaxIntSkew;  // Max Internal Skewness
  /** @brief Max Concativity
   **/
  double qcMaxConc;  // Max Concativity
  /** @brief Minimum Cell Volume
   **/
  double qcMinVol;  // Minimum Cell Volume
  /** @brief Minimum Tet Quality
   **/
  double qcMinTetQ;  // Minimum Tet Quality
  /** @brief Minimum Area
   **/
  double qcMinArea;  // Minimum Area
  /** @brief Minimum Twist
   **/
  double qcMinTwist;  // Minimum Twist
  /** @brief Minimum Face Weight
   **/
  double qcMinFaceW;  // Minimum Face Weight
  /** @brief Minimum Volume Ratio
   **/
  double qcMinVolRto;  // Minimum Volume Ratio
  /** @brief Minimum Determinant
   **/
  double qcMinDet;  // Minimum Determinant
  /** @brief Minimum Triangle Twist
   **/
  double qcMinTrTwist;  // Minimum Triangle Twist
  /** @brief nSmoothScale
   **/
  int qcSmthScale;  // nSmoothScale
  /** @brief Error Reduction
   **/
  double qcErrRedctn;  // Error Reduction

  // --- Misc. General
  /** @brief merge tolerance for mesh
   **/
  double mergeTol;  // Merge Tolerance

  /** @brief Vector for shmRegionRefine struct
   **/
  std::vector<shmRegionRefine> geomRefs;

  /** @brief Vector for shmSurfRefine struct
   **/
  std::vector<shmSurfRefine> surfRefs;

  /** @brief Vector for STL patches
   **/
  std::vector<shmSTLDefinition> stlPatchDefs;

  /** @brief Vector for searchable shapes
   **/
  std::vector<shmSearchableShapes> srchShape;

  /** @brief Vector for feature edge refinement
   **/
  std::vector<shmFeatureEdgeRef> ftrEdge;

  /** @brief Vector for Layers
   **/
  std::vector<shmLayers> layerVec;
};

#endif

#endif
