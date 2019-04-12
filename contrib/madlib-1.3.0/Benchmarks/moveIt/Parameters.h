// -*- C++ -*-
// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
//
// Authors: Gaetan Compere, Jean-Francois Remacle
// -------------------------------------------------------------------

#ifndef _H_PARAMETERS
#define _H_PARAMETERS

#include "AdaptInterface.h"
#include "SizeFieldBase.h"
#include "LocalSizeField.h"

#include <string.h>
using std::string;
#include <list>
#include <set>
#include <vector>
using std::vector;

// ----------------------------------------------------------------------
enum KinematicType {
  KT_DISPLACEMENT,
  KT_VELOCITY
};

// ----------------------------------------------------------------------
enum OrientationType {
  ORT_ISOTROPIC,
  ORT_ANISOTROPIC
};

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
struct LocalSFDef {
  string name;
  bool distToFaces;
  bool isotropic;
  double radius;
  string sizeN, sizeT;
  bool limit;
  double limitSizeTg;
  double maxCurv;
};

// ----------------------------------------------------------------------
struct SizeFieldDef
{
  MAd::sFieldType type; // analytical or piecewise linear

  // analytical size fields
  OrientationType orientation;
  string isoSize;
  vector<string> anisoSize;
  vector<string> anisoDir0, anisoDir1, anisoDir2;

  // piecewise linear size fields
  string pwlSource; // curvature, initial length or file
  
  // curvature parameters
  bool curv_aniso;
  double curv_alpha, curv_hMin;

  // file parameters
  string mshFile;
};

// ----------------------------------------------------------------------
struct ObjectDef {

  string name;

  // geometry
  int geomLevel; // 0=point,1=line,2=surface,3=region
  std::set<int> geomTags;

  // kinematics
  KinematicType kinType;
  vector<string> kinExpression;

  // size fields
  std::list<LocalSFDef> sizes;
};

// ----------------------------------------------------------------------
struct MoveParameters {
  
  // global input specifications
  // ---------------------------

  string meshFileName;
  string geoFileName;
  
  // global output parameters
  // ------------------------

  int outFrequency;
  string outType; // "Msh", "Pos" or "MshAndPos"(def)
  string outPrefix;

  // Global task parameters
  // ----------------------

  string task; // MobileObject(def),...

  // debugging parmeters
  // --------------------

  int debugLevel;
  bool openJournal;
  string referenceJournal;
  bool sliverReports;
  bool testSliverOperators;

  // Mesh adaptation
  // ----------------

  bool adapt;
  bool adaptAlways;
  double adaptQualityThr;
  int maxInnerIter;
  double lowerLength, upperLength;
  double infLength;
  double swapMinImproveRatio;
  double sliverQuality;

  bool splitCanMakeSlivers, collapseCanMakeSlivers;
  double makeSliverInSplitLenSqBound, makeSliverInCollapseLenSqBound;

  bool collapseOnBoundary;
  double clpOnBdryTolerance;

  bool swapOnBoundary;
  double swapOnBdryTolerance;

  bool trackGeometry;
  bool   snap_cavityIsMesh;
  int    snap_thickness;
  double snap_chi;
  bool   snap_check;
  bool   snap_force;

  bool SFSmoothing;
  double SFSmoothGrad;

  // Moving objects
  // --------------

  std::vector<ObjectDef> objects;

  // Size Field parameters
  // ---------------------

  std::list<SizeFieldDef> sizes;
  
  // Node motion
  // -----------

  string nodeMotionType; // None(def), ForceBoundaries, ElasticAnalogy
  bool   elasticSubAdapt;
  double elasticSubAdaptQualThr;
  bool   elasticIsMeshTheCavity;
  int    elasticCavityThickness;
  double elasticStiffnessAlteration;

  // Time parameters
  // ---------------

  int maxNbTimeSteps;
  double timeStep;
  double finalTime;

};

// ----------------------------------------------------------------------
// this function set the parameters to the values of a particular test case
extern void setCurrentParameters(MoveParameters *);

// ----------------------------------------------------------------------
void setDefaultParameters(MoveParameters * params)
{
  params->meshFileName = "";
  params->geoFileName = "";

  params->outFrequency = 1;
  params->outType = "MshAndPos";
  params->outPrefix  = "result/";

  params->task = "MobileObject";

  params->debugLevel = 1;
  params->openJournal = false;
  params->referenceJournal = "";
  params->sliverReports = false;
  params->testSliverOperators = false;

  params->adapt = true;
  params->adaptAlways = true;
  params->adaptQualityThr = 0.02;
  
  params->maxInnerIter = 10;
  params->lowerLength = 1./sqrt(3.);
  params->upperLength = sqrt(3.);
  params->infLength = 1.e14;
  params->swapMinImproveRatio = 1.0;
  params->sliverQuality = 0.02;

  params->splitCanMakeSlivers = true;
  params->makeSliverInSplitLenSqBound = 10.;
  params->collapseCanMakeSlivers = true;
  params->makeSliverInCollapseLenSqBound = 0.1;

  params->collapseOnBoundary = true;
  params->clpOnBdryTolerance = 1.e-6;

  params->swapOnBoundary = true;
  params->swapOnBdryTolerance = 1.e-6;

  params->trackGeometry = false;
  params->snap_cavityIsMesh = false;
  params->snap_thickness = 3;
  params->snap_chi = 1.;
  params->snap_check = false;
  params->snap_force = false;

  params->SFSmoothing = false;
  params->SFSmoothGrad = 1.;

  params->objects.clear();

  params->sizes.clear();

  params->nodeMotionType = "None";
  params->elasticSubAdapt = true;
  params->elasticSubAdaptQualThr = 0.;
  params->elasticIsMeshTheCavity = true;
  params->elasticCavityThickness = 3;
  params->elasticStiffnessAlteration = 1.0;

  params->maxNbTimeSteps = 1000000;
  params->timeStep = 0.1;
  params->finalTime = -1.0;
}

// ----------------------------------------------------------------------
MoveParameters getParameters()
{
  MoveParameters params;
  setDefaultParameters(&params);
  setCurrentParameters(&params);

  return params;
}

// ----------------------------------------------------------------------

#endif
