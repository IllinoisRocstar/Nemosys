#ifndef SNAPPYMESHPARAMS_H
#define SNAPPYMESHPARAMS_H

#ifdef HAVE_CFMSH 

#include "meshingParams.H"
#include <string>
#include <vector>
#include <map>
#include <tuple>



struct shmGeomRefine // Define Geometry
{
  std::string patchNm;
  std::string searchableName;
  std::string shapeParameters1; // Points, Center etc
  std::string shapeParameters2; // Points, Center etc
  double rad;           // Radius
  std::string mode;       // Inside, Outside
  int minLvl;
  int maxLvl;

};


struct shmRegionRef // Region Refinement (Add all shapes and their controls)
{
  std::string refPatchNm;
  int minLvl;
  int maxLvl;
};

/*struct shmLayers // Layer Addition on Patches
{

}*/ 


class snappymeshParams : public meshingParams
{
  public:
  // initializing params with defaults
  snappymeshParams();
  ~snappymeshParams(){};


  // Booleans
  bool _withCastMesh;     // Boolean for castellated mesh
  bool _withSnap;       // Boolean for snap
  bool _withLayers;     // Boolean for add layers
  bool _withGeomRefReg;       // Geometry refinement regions (searchable)
  bool _withCellZones;      // Defining cellzones and facezones enabled?
  bool _withSurfRefReg;   // Surface refinement regions

  // Castellated Mesh Controls
  std::string geomFileName; // Input geometry STL name
  std::string geomPatchName;  // patch name for same geometry file
  std::string surfRefPatch; // surface refinement patch stl
  int maxLCells;        // max global cells
  int maxGCells;        // max local cells (on 1 processor)
  int minRefCells;      // minimum refinement cells
  int cellsBetnLvls;      // number of cells between levels
  int refSurfLvlMin;      // minimum surface refinement
  int refSurfLvlMax;      // maximum surface refinement
  double featAngle;     // resolve feature angle
  double locMeshX;      // location in Mesh
  double locMeshY;      // location in Mesh
  double locMeshZ;      // location in Mesh
  bool _alwFreeZone;      // allow free standing zones


  // Snap Controls
  int snapSmthPatch;      // nSmoothPatch
  double snapTol;       // Snap Tolerance
  int solveSnapIter;      // nSolveIter
  int relaxSnapIter;      // nRelaxIter


  // Layer Controls
  bool _relSize;        // Relative Sizes
  double expRatio;      // Expansion Ratio
  double finLThick;     // Final Layer Thickness
  double minThick;      // Minimum Thickness
  int nGrow;
  double lyrFeatAngle;    // Feature Angle
  int lyrRelaxIter;     // number of relaxation interations
  int lyrSmthSurfNorm;    // # of smooth surface normals
  int lyrSmthNorm;      // # of smooth normals
  int lyrSmthThick;     // # of smooth thickness
  double lyrMaxFcTR;      // Maximum face thickness ratio
  double lyrMaxThickTMR;    // Maximum thickness to medial ratio
  double lyrMinMedAngl;     // Minimum medial axis angle
  int lyrBuffrCells;      // # of buffer cells no extrude
  int lyrIter;        // # of layer interations


  // Mesh Quality Controls
  int qcMaxNOrtho;      // Maximum non-orthogonality
  double qcMaxBndrySkew;      // Max Boundary Skewness
  double qcMaxIntSkew;      // Max Internal Skewness
  double qcMaxConc;       // Max Concativity
  double qcMinVol;      // Minimum Cell Volume
  double qcMinTetQ;     // Minimum Tet Quality
  double qcMinArea;     // Minimum Area
  double qcMinTwist;      // Minimum Twist
  double qcMinFaceW;      // Minimum Face Weight
  double qcMinVolRto;     // Minimum Volume Ratio
  double qcMinDet;      // Minimum Determinant
  double qcMinTrTwist;    // Minimum Triangle Twist
  int qcSmthScale;      // nSmoothScale
  double qcErrRedctn;     // Error Reduction


  // Misc. General
  double mergeTol;      // Merge Tolerance

  std::vector<shmGeomRefine> geomRefs;
  std::vector<shmRegionRef> surfRefs;

};

#endif

#endif