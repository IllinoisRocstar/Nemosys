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

#ifndef _H_MYPARAMS
#define _H_MYPARAMS

#include "Parameters.h"

// ----------------------------------------------------------------------
void setCurrentParameters(MoveParameters * params) {

  params->meshFileName = "tube.msh";
//   params->geoFileName="tube_fullGeo.geo";
  
  // global output parameters
  // ------------------------

  params->outFrequency = 10;
  params->outType = "Msh"; // "Msh", "Pos" or "MshAndPos"
  params->outPrefix = "result/";

  // Global task parameters
  // ----------------------

  params->task = "MobileObject";

  // debugging parmeters
  // --------------------

  params->debugLevel = 0;
  params->openJournal = false;
  params->referenceJournal = "";
  params->sliverReports = false;
  params->testSliverOperators = false;

  // Mesh adaptation
  // ----------------

//   params->maxInnerIter = 10;
//   params->lowerLength = 0.577;
//   params->upperLength = 1.732;
//   params->swapMinImproveRatio = 1.0;
//   params->sliverQuality = 0.02;
//   params->splitCanMakeSlivers = true;
//   params->makeSliverInSplitLenSqBound = 0.1;
//   params->collapseCanMakeSlivers = true;
//   params->makeSliverInCollapseLenSqBound = 10.;
//   params->collapseOnBoundary = true;
//   params->clpOnBdryTolerance = 1.e-6;
//   params->swapOnBoundary = true;
//   params->swapOnBdryTolerance = 1.e2;
//   params->trackGeometry = false;
//   params->snap_cavityIsMesh = false;
//   params->snap_thickness = 3;
//   params->snap_chi = 1.;
//   params->snap_check = false;
//   params->snap_force = false;
  
  params->SFSmoothing = true;
  params->SFSmoothGrad = 1.;

  // Node motion
  // -----------

  params->nodeMotionType = "ElasticAnalogy";
  params->elasticSubAdapt = true;
  params->elasticSubAdaptQualThr = 0.;
  params->elasticStiffnessAlteration = 1.0;
  params->elasticIsMeshTheCavity = true;
  params->elasticCavityThickness = 0;

  // Time parameters
  // ---------------

  params->timeStep = 0.01;
  params->finalTime = 4.0;
//   params->maxNbTimeSteps = 1;

  
  // Moving objects
  // --------------

  // --- The cylinder ---

  ObjectDef cylinder;
  cylinder.name = "Cylinder";
  cylinder.geomLevel = 2;
  cylinder.geomTags.insert(8);
  cylinder.geomTags.insert(12);
  cylinder.geomTags.insert(16);
  cylinder.geomTags.insert(20);
  cylinder.geomTags.insert(22);
  cylinder.geomTags.insert(126);
  cylinder.kinType = KT_VELOCITY;
  cylinder.kinExpression.push_back("0.");
  cylinder.kinExpression.push_back("0.");
  cylinder.kinExpression.push_back("1.*sign(sin(2.*3.14159265/8.*t+1.e-6))");

  LocalSFDef localSFCyl;
  localSFCyl.name = "CylinderSF";
	localSFCyl.radius = 1.;
  localSFCyl.distToFaces = false;
  localSFCyl.isotropic = true;
  localSFCyl.sizeN = "0.1+0.5*x";

  cylinder.sizes.push_back(localSFCyl);

  params->objects.push_back(cylinder);

  // --- The tube ---

  ObjectDef tube;
  tube.name = "Tube";
  tube.geomLevel = 2;
  tube.geomTags.insert(40);
  tube.geomTags.insert(44);
  tube.geomTags.insert(48);
  tube.geomTags.insert(52);
  tube.geomTags.insert(56);
  tube.geomTags.insert(60);
  tube.geomTags.insert(64);
  tube.geomTags.insert(68);
  tube.geomTags.insert(71);
  tube.geomTags.insert(106);
  tube.kinType = KT_VELOCITY;
  tube.kinExpression.push_back("0.");
  tube.kinExpression.push_back("0.");
  tube.kinExpression.push_back("-1.*sign(sin(2.*3.14159265/8.*t+1.e-6))");

  LocalSFDef localSFTube;
  localSFTube.name = "tubeSF";
	localSFTube.radius = 1.;
  localSFTube.distToFaces = false;
  localSFTube.isotropic = true;
  localSFTube.sizeN = "0.2+0.4*x";

  tube.sizes.push_back(localSFTube);

  params->objects.push_back(tube);

  // --- The box ---

//   ObjectDef box;
//   box.name = "Box";
//   box.geomLevel = 2;
//   box.geomTags.insert(110);
//   box.geomTags.insert(114);
//   box.geomTags.insert(118);
//   box.geomTags.insert(122);
//   box.geomTags.insert(124);
//   box.geomTags.insert(131);
//   box.kinType = KT_DISPLACEMENT;
//   box.kinExpression.push_back("0.");
//   box.kinExpression.push_back("0.");
//   box.kinExpression.push_back("0.");

//   params->objects.push_back(box);

  // Size Field parameters
  // ---------------------

  SizeFieldDef analSF;
  analSF.type = MAd::ANALYTICALSFIELD;
  analSF.orientation = ORT_ISOTROPIC;
  analSF.isoSize = "0.6";
//   analSF.orientation = ORT_ANISOTROPIC;
//   vector<string> myAnisoSize;
//   myAnisoSize.push_back("0.3");
//   myAnisoSize.push_back("0.6");
//   myAnisoSize.push_back("0.6");
//   analSF.anisoSize = myAnisoSize;
//   vector<string> anisoDir0;
//   anisoDir0.push_back("1.");
//   anisoDir0.push_back("0.");
//   anisoDir0.push_back("0.");
//   analSF.anisoDir0 = anisoDir0;
//   vector<string> anisoDir1;
//   anisoDir1.push_back("0.");
//   anisoDir1.push_back("1.");
//   anisoDir1.push_back("0.");
//   analSF.anisoDir1 = anisoDir1;
//   vector<string> anisoDir2;
//   anisoDir2.push_back("0.");
//   anisoDir2.push_back("0.");
//   anisoDir2.push_back("1.");
//   analSF.anisoDir2 = anisoDir2;

  params->sizes.push_back(analSF);
}

// ----------------------------------------------------------------------

#endif
