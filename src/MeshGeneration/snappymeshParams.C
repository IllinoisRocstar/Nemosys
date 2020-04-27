#include "snappymeshParams.H"

// Add default snappyHexMesh parameters here

snappymeshParams::snappymeshParams() {
  // All default values are defined here for initialization
  /*_withCastMesh = true;     // Boolean for castellated mesh
  _withSnap = true;       // Boolean for snap
  _withLayers = true;     // Boolean for add layers
  _withCellZones = false;*/

  // Geometry definition
  _withMultiPatches = false;
  geomTolerance = 1e-05;  // Keep optional
  maxTreeDepth = 10;      // Keep optional

  // // Castellated Mesh Controls
  maxLCells = 2000000;  // max global cells
  maxGCells = 4000000;  // max local cells (on 1 processor)
  minRefCells = 0;      // minimum refinement level
  cellsBetnLvls = 3;    // number of cells between levels
  refSurfLvlMin = 0;    // minimum surface refinement
  refSurfLvlMax = 0;    // maximum surface refinement
  featAngle = 60;       // resolve feature angle
  locMeshX = 0;         // location in Mesh
  locMeshY = 0;         // location in Mesh
  locMeshZ = 0;         // location in Mesh
  _alwFreeZone = true;  // allow free standing zones
  gPLvlInc = 1;         // Gap level increment
  planarAngle = 30;
  castMeshGpLvl = 1;

  // Snap Controls
  snapSmthPatch = 4;               // nSmoothPatch
  snapTol = 0.5;                   // Snap Tolerance
  solveSnapIter = 200;             // nSolveIter
  relaxSnapIter = 6;               // nRelaxIter
  nFeatureSnapIter = 10;           // nFeatureSnapIter
  implicitFeatureSnap = false;     // implicitFeatureSnap
  explicitFeatureSnap = false;     // explicitFeatureSnap
  multiRegionFeatureSnap = false;  // multiRegionFeatureSnap

  // Layer Controls
  _relSize = true;        // Relative Sizes
  expRatio = 1.3;         // Expansion Ratio
  finLThick = 1;          // Final Layer Thickness
  minThick = 0.1;         // Minimum Thickness
  nGrow = 0;              // nGrow
  lyrFeatAngle = 30;      // Feature Angle
  lyrRelaxIter = 3;       // number of relaxation interations
  lyrSmthSurfNorm = 1;    // # of smooth surface normals
  lyrSmthNorm = 3;        // # of smooth normals
  lyrSmthThick = 2;       // # of smooth thickness
  lyrMaxFcTR = 0.5;       // Maximum face thickness ratio
  lyrMaxThickTMR = 1;     // Maximum thickness to medial ratio
  lyrMinMedAngl = 90;     // Minimum medial axis angle
  lyrBuffrCells = 0;      // # of buffer cells no extrude
  lyrIter = 50;           // # of layer interations
  nRelaxedIter = 20;      // # of relaxed iteration
  slipFeatureAngle = 30;  // slipFeatureAngle

  // Mesh Quality Controls
  qcMaxNOrtho = 65;     // Maximum non-orthogonality
  qcMaxBndrySkew = 20;  // Max Boundary Skewness
  qcMaxIntSkew = 4;     // Max Internal Skewness
  qcMaxConc = 80;       // Max Concativity
  qcMinVol = 1e-13;     // Minimum Cell Volume
  qcMinTetQ = 1e-15;    // Minimum Tet Quality
  qcMinArea = -1;       // Minimum Area
  qcMinTwist = 0.02;    // Minimum Twist
  qcMinFaceW = 0.05;    // Minimum Face Weight
  qcMinVolRto = 0.01;   // Minimum Volume Ratio
  qcMinDet = 0.001;     // Minimum Determinant
  qcMinTrTwist = -1;    // Minimum Triangle Twist
  qcSmthScale = 5;      // nSmoothScale
  qcErrRedctn = 0.75;   // Error Reduction

  // Misc. General
  mergeTol = 1e-6;  // Merge Tolerance
}
