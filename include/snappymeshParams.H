#ifndef SNAPPYMESHPARAMS_H
#define SNAPPYMESHPARAMS_H

#ifdef HAVE_CFMSH 

#include "meshingParams.H"
#include <string>
#include <vector>
#include <map>
#include <tuple>


/**
    @brief A structure to respresent geometry patches
          refining capability in snappymeshGen.
**/
struct shmGeomRefine // Define Geometry
{
  std::string patchNm;  /**< Patch name */
  std::string searchableName; /**< Patch name */
  std::string shapeParameters1; /**< Points, Center etc */
  std::string shapeParameters2; /**< Points, Center etc */
  double rad;           /**< Radius */
  std::string mode;       /**< Inside, Outside */
  int minLvl; /**< Minimum refinement level */
  int maxLvl; /**< Maximum refinement level */

};

/**
    @brief A structure to represent geometry regions
          refining capability in snappymeshgGen
**/
struct shmRegionRef // Region Refinement (Add all shapes and their controls)
{
  std::string refPatchNm; /**< Refining region patch name */
  int minLvl; /**< Minimum refinement level */
  int maxLvl; /**< Maximum refinement level */
};

// Future development module
/*struct shmLayers // Layer Addition on Patches
{

}*/ 



/** 
  @brief snappymeshParams contains all parameters essential for mesh generation
      using snappymeshGen class methods. These parameters are assigned a value
      at time of user input parsing through JSON in meshGemDriver.
**/
class snappymeshParams : public meshingParams
{
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
  bool _withCastMesh;     // Boolean for castellated mesh
  
  /** @brief Enables surface snapping option for snappymeshGen
  **/
  bool _withSnap;       // Boolean for snap
  
  /** @brief Enables add layers option for snappymeshGen
  **/
  bool _withLayers;     // Boolean for add layers
  
  /** @brief Enables user defined geometry refinement regions (searchable)
  **/
  bool _withGeomRefReg;       // Geometry refinement regions (searchable)
  
  /** @brief Enables defining cellzones and facezones
  **/
  bool _withCellZones;      // Defining cellzones and facezones enabled?
  
  /** @brief Enables user defined surface refinement regions
  **/
  bool _withSurfRefReg;   // Surface refinement regions


  // --- Castellated Mesh Controls
  /** @brief Input geometry STL name
  **/
  std::string geomFileName; // Input geometry STL name
  
  /** @brief patch name for input surface file
  **/
  std::string geomPatchName;  // patch name for same geometry file
  
  /** @brief surface refinement patch surface file
  **/
  std::string surfRefPatch; // surface refinement patch stl
  
  /** @brief max global cells allowed in mesh
  **/
  int maxLCells;        // max global cells
  
  /** @brief max local cells allowed (on 1 processor)
  **/
  int maxGCells;        // max local cells (on 1 processor)
  
  /** @brief minimum refinement cells
  **/
  int minRefCells;      // minimum refinement cells
  
  /** @brief number of cells between levels
  **/
  int cellsBetnLvls;      // number of cells between levels
  
  /** @brief minimum surface refinement
  **/
  int refSurfLvlMin;      // minimum surface refinement
  
  /** @brief maximum surface refinement
  **/
  int refSurfLvlMax;      // maximum surface refinement
  
  /** @brief resolve feature angle
  **/
  double featAngle;     // resolve feature angle
  
  /** @brief location in Mesh (cells are kept in regions accessinble from this
            location)
  **/
  double locMeshX;      // location in Mesh
  
  /** @brief location in Mesh (cells are kept in regions accessinble from this
            location)
  **/
  double locMeshY;      // location in Mesh
  
  /** @brief location in Mesh (cells are kept in regions accessinble from this
            location)
  **/
  double locMeshZ;      // location in Mesh
  
  /** @brief allows free standing zones if enabled
  **/
  bool _alwFreeZone;      // allow free standing zones


  // --- Snap Controls
  /** @brief Number of smoothing patches during snapping procedure
  **/
  int snapSmthPatch;      // nSmoothPatch
  
  /** @brief Snapping tolerance at shared interface (affects conformality of 
            mesh)
  **/
  double snapTol;       // Snap Tolerance
  
  /** @brief Maximum iterations during snapping procedure before ending the 
            process
  **/
  int solveSnapIter;      // nSolveIter
  
  /** @brief Maximum number of relaxation iterations for snapping procedure
  **/
  int relaxSnapIter;      // nRelaxIter


  // --- Layer Controls
  /** @brief Enables relative sizes option during layer addition
  **/
  bool _relSize;        // Relative Sizes
  
  /** @brief Expansion Ratio for layer addition
  **/
  double expRatio;      // Expansion Ratio
  
  /** @brief Final Layer Thickness for layer addition
  **/
  double finLThick;     // Final Layer Thickness
  
  /** @brief Minimum Thickness for layer addition
  **/
  double minThick;      // Minimum Thickness
  
  /** @brief Growth rate of successive layers
  **/
  int nGrow;
  
  /** @brief Feature Angle for layer addition
  **/
  double lyrFeatAngle;    // Feature Angle
  
  /** @brief number of relaxation interations for layer addition
  **/
  int lyrRelaxIter;     // number of relaxation interations
  
  /** @brief Number of smooth surface normals for layer addition
  **/
  int lyrSmthSurfNorm;    // # of smooth surface normals
  
  /** @brief Number of smooth normals for layer addition
  **/
  int lyrSmthNorm;      // # of smooth normals
  
  /** @brief Number of smooth thickness for layer addition
  **/
  int lyrSmthThick;     // # of smooth thickness
  
  /** @brief Maximum face thickness ratio for layer addition
  **/
  double lyrMaxFcTR;      // Maximum face thickness ratio
  
  /** @brief Maximum thickness to medial ratio for layer addition
  **/
  double lyrMaxThickTMR;    // Maximum thickness to medial ratio
  
  /** @brief Minimum medial axis angle for layer addition
  **/
  double lyrMinMedAngl;     // Minimum medial axis angle
  
  /** @brief Number of buffer cells no extrude for layer addition
  **/
  int lyrBuffrCells;      // # of buffer cells no extrude
  
  /** @brief Number of layer interations for layer addition
  **/
  int lyrIter;        // # of layer interations


  // --- Mesh Quality Controls
  /** @brief Maximum non-orthogonality
  **/
  int qcMaxNOrtho;      // Maximum non-orthogonality
  /** @brief Max Boundary Skewness
  **/
  double qcMaxBndrySkew;      // Max Boundary Skewness
  /** @brief Max Internal Skewness
  **/
  double qcMaxIntSkew;      // Max Internal Skewness
  /** @brief Max Concativity
  **/
  double qcMaxConc;       // Max Concativity
  /** @brief Minimum Cell Volume
  **/
  double qcMinVol;      // Minimum Cell Volume
  /** @brief Minimum Tet Quality
  **/
  double qcMinTetQ;     // Minimum Tet Quality
  /** @brief Minimum Area
  **/
  double qcMinArea;     // Minimum Area
  /** @brief Minimum Twist
  **/
  double qcMinTwist;      // Minimum Twist
  /** @brief Minimum Face Weight
  **/
  double qcMinFaceW;      // Minimum Face Weight
  /** @brief Minimum Volume Ratio
  **/
  double qcMinVolRto;     // Minimum Volume Ratio
  /** @brief Minimum Determinant
  **/
  double qcMinDet;      // Minimum Determinant
  /** @brief Minimum Triangle Twist
  **/
  double qcMinTrTwist;    // Minimum Triangle Twist
  /** @brief nSmoothScale
  **/
  int qcSmthScale;      // nSmoothScale
  /** @brief Error Reduction
  **/
  double qcErrRedctn;     // Error Reduction


  // --- Misc. General
  /** @brief merge tolerance for mesh 
  **/
  double mergeTol;      // Merge Tolerance

  /** @brief Vector for shmGeomRefine struct
  **/
  std::vector<shmGeomRefine> geomRefs;

  /** @brief Vector for shmRegionRef struct
  **/
  std::vector<shmRegionRef> surfRefs;

};

#endif

#endif