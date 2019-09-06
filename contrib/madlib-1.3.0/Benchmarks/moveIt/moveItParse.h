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

#ifndef _H_MOVEITPARSE_
#define _H_MOVEITPARSE_

#ifdef _HAVE_PARSER_

#include "Parser.h"

#include "MAdLib.h"
#include "Parameters.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

extern int yyCMDparse();

//==============================================================================
// Control file definition 
//==============================================================================

void DefineMMCtrlFile() {
  
  // global input specifications
  // ---------------------------

  ClassParameter& input = Parser::instance().add("Input");
  ClassParameter& inputMesh = input.add("Mesh");
  inputMesh.addString("MshFile","");
  inputMesh.addString("GeoFile","");
  
  // global output parameters
  // ------------------------
  
  ClassParameter& output = Parser::instance().add("Output");
  output.addString("Prefix","result/");
  output.addInteger("Frequency",1);
  output.addToken("Type",3,"MshAndPos","Msh","Pos","MshAndPos");

  // Task parameters
  // ----------------

  ClassParameter& task = Parser::instance().add("Task");
  task.addToken("Type",24,
                "MobileObject",
                "RemoveRegion",
                "SizeField",
                "RemoveSliverFaces",
                "LaplaceSmoothing",
                "GeoConstraints",
                "TopoConstraints",
                "FaceSwap",
                "UglyMesh",
                "RemoveSliverRegions",
                "DESC",
                "EdgeSwap",
                "GeoMatcher",
                "EdgeSplit",
                "EdgeCollapse",
                "EdgeSplitCollapse",
                "EdgeSplitCollapse2",
                "FaceCollapse",
                "OptimizeLength",
                "OptimizeShape",
                "Metric",
                "Curvatures",
                "Miscellaneous",
                "None",
                "MobileObject");

  // debug parameters
  // -----------------
  
  ClassParameter& debug = Parser::instance().add("Debug");
  debug.addInteger("DebugLevel",1);
  debug.addToken("OpenJournal",2,"Yes","No","No");
  debug.addString("ReferenceJournal","");
  debug.addToken("ReportSlivers",2,"Yes","No","No");
  debug.addToken("TestSliverOperators",2,"Yes","No","No");

  // Mesh adaptation
  // ----------------

  ClassParameter& ma = Parser::instance().add("Adaptation");
  ma.addToken("Enable",2,"No","Yes","Yes");
  ma.addToken("AdaptEveryIteration",2,"No","Yes","Yes");
  ma.addDouble("QualityThreshold",0.02);
  ma.addInteger("MaxNumberOfInnerIteration",10);
  ma.addDouble("LowerLengthBound",1./sqrt(3.));
  ma.addDouble("UpperLengthBound",sqrt(3.));
  ma.addDouble("InfiniteEdgeLength",1.e14);
  ma.addDouble("SwapMinImproveRatio",1.0);
  ma.addDouble("SliverQuality",0.02);
  ma.addToken("SplitCanMakeSlivers",2,"No","Yes","Yes");
  ma.addToken("CollapseCanMakeSlivers",2,"No","Yes","Yes");
  ma.addDouble("SliverLowerLengthSqBound",0.1);
  ma.addDouble("SliverUpperLengthSqBound",10.);
  ma.addToken("CollapseOnBoundary",2,"No","Yes","Yes");
  ClassParameter& mabound = ma.add("BoundaryCollapses");
  mabound.addDouble("Tolerance",1.e-6);
  ma.addToken("SwapOnBoundary",2,"No","Yes","Yes");
  ClassParameter& maboundswp = ma.add("BoundarySwaps");
  maboundswp.addDouble("Tolerance",1.e-6);
  ma.addToken("TrackGeometry",2,"No","Yes","No");
  ClassParameter& masnap = ma.add("VertexSnapping");
  masnap.addToken("ForceRelocation",2,"No","Yes","No");
  masnap.addToken("StrictChecking",2,"No","Yes","No");
  masnap.addToken("IsMeshTheCavity",2,"No","Yes","No");
  masnap.addInteger("CavityThickness",3);
  masnap.addDouble("StiffnessAlteration",1.0);
  ma.addToken("SmoothSizeField",2,"No","Yes","No");
  ClassParameter& smoo = ma.add("SizeFieldSmoothing");
  smoo.addDouble("MaximumGradient",1.0);
  
  // moving objects
  // --------------

  ClassParameter& mo = Parser::instance().add("MovingObjects");
  ClassParameter& object = mo.addMultiple("Object");

  // geometry
  ClassParameter& obj_geom = object.add("Geometry");
  obj_geom.addToken("GeometricalType",4,"Vertex","Edge","Face","Region","Face");
  obj_geom.addIntegerRange("GeometricalTags");

  // kinematics
  ClassParameter& obj_kin = object.add("Kinematics");
  obj_kin.addToken("KinematicsSpecification",2,"Displacement","Velocity","Displacement");
  ClassParameter& kinDx = obj_kin.add("Displacement");
  kinDx.addStringList("Formulation");
  ClassParameter& kinV = obj_kin.add("Velocity");
  kinV.addStringList("Formulation");

  // size field
  ClassParameter& locS = object.addMultiple("LocalSizeField");
//   locS.addDouble("Radius",0.);
//   locS.addToken("SizeFunction",2,"Linear","Sqrt","Linear");
//   locS.addDouble("SizeOnObject",1.);
//   locS.addDouble("SizeOutside",1.);
//   locS.addToken("Orientation",2,"Isotropic","Anisotropic","Isotropic");
//   ClassParameter& anisoLsf = locS.add("Anisotropic");
//   anisoLsf.addDouble("SizeTangent",1.);
  locS.addToken("DistanceComputation",2,"ToWallVertices","ToWallFaces","ToWallVertices");
  locS.addToken("Orientation",2,"Isotropic","Anisotropic","Isotropic");
  locS.addDouble("Radius",0.);
  ClassParameter& isoLsf = locS.add("Isotropic");
  isoLsf.addString("Size","");
  ClassParameter& anisoLsf = locS.add("Anisotropic");
  anisoLsf.addString("NormalSize","");
  anisoLsf.addString("TangentSize","");
  ClassParameter& limiter = locS.add("CurvatureBasedLimiter");
  limiter.addToken("Enable",2,"Yes","No","No");
  limiter.addDouble("TangentSizeOverRadius",1.);
  limiter.addDouble("CurvatureLimiter",1.e6);

  // Size Field parameters
  // ---------------------

  ClassParameter& sf = Parser::instance().addMultiple("SizeField");
  sf.addToken("Type",3,"Analytical","PiecewiseLinear","Background","Analytical");
  ClassParameter& sf_an = sf.add("Analytical");
  sf_an.addToken("OrientationType",2,"Isotropic","Anisotropic","Isotropic");
  ClassParameter& sf_an_iso = sf_an.add("Isotropic");
  sf_an_iso.addString("Length","");
  ClassParameter& sf_an_aniso = sf_an.add("Anisotropic");
  sf_an_aniso.addStringList("Lengths");
  sf_an_aniso.addStringList("Direction1");
  sf_an_aniso.addStringList("Direction2");
  sf_an_aniso.addStringList("Direction3");
  ClassParameter& sf_pwl = sf.add("PiecewiseLinear");
  sf_pwl.addToken("Source",2,"InitialLength","Curvature","Curvature");
  ClassParameter& sf_pwl_curv = sf_pwl.add("Curvature");
  sf_pwl_curv.addToken("Anisotropic",2,"Yes","No","Yes");
  sf_pwl_curv.addDouble("Alpha",2.);
  sf_pwl_curv.addDouble("MinimalLength",1.e-4);
  ClassParameter& sf_bg = sf.add("Background");
  sf_bg.addString("FileName","");
  
  // Node motion
  // -----------

  ClassParameter& nm = Parser::instance().add("NodeMotion");
  nm.addToken("Type",3,"ElasticAnalogy","ForceBoundaries","None","None");
  ClassParameter& maelas = nm.add("ElasticParameters");
  maelas.addToken("AllowSubAdaptation",2,"No","Yes","Yes");
  maelas.addDouble("SubAdaptationQualityThreshold",0.0);
  maelas.addToken("IsMeshTheCavity",2,"No","Yes","Yes");
  maelas.addInteger("CavityThickness",3);
  maelas.addDouble("StiffnessAlteration",1.0);

  // Time parameters
  // ---------------

  ClassParameter& time = Parser::instance().add("TimeSpecifications"); 
  time.addInteger("MaxTimeSteps",1000000);
  time.addDouble("TimeStep",0.1);
  time.addDouble("FinalTime",-1.0);
}

//==============================================================================
// Parsing functions
//==============================================================================

string parseMMCommandLine(int argc,char** argv) {

  string inputFile = "";
  string templateFile = "";

  bool inputFileSpecified = false;
  bool dumpTemplate   = false;

  if (argc == 1) {
    dumpTemplate = true;
  }
  else if ( argc==2 ){
    inputFileSpecified = true;
    inputFile = argv[1];
  }
  else {
    std::cerr << "Usages:" << argv[0]<<"  control_file"<<std::endl;
    exit(1);
  }
  
  // dump a template file 
  
  if (dumpTemplate) {
    templateFile = "template.mad";
    FILE* fpcheck = fopen(templateFile.c_str(),"w");
    fprintf(fpcheck,"//-*- C++ -*-\n");
    printf("Writing a template parameter tree to file \'%s\'\n",templateFile.c_str());
    Parser::instance().print(fpcheck);  
    fclose(fpcheck);
    exit(1);
  } 

  // read the parameters 
  
  if (!inputFileSpecified || inputFile == "") {
    printf("No parameter file was specified - exiting");
    exit(1);
  }
  
  return inputFile;
}

// -----------------------------------------------------------------------------
bool parseMMCtrlFile(const string& inputFile) {

  FILE* fp = freopen(inputFile.c_str(), "r", stdin);
  if (!fp) {
    printf("Error: could not open parameter file \'%s\'\n",inputFile.c_str());
    return false;
  }
  
  yyCMDparse();

  return true;
}

// ----------------------------------------------------------------------
// Parse the control file
void parseControlFile(int argc, char* argv[])
{
  cout << "Parsing the control file...\n";
  DefineMMCtrlFile();
  string inputFile = parseMMCommandLine(argc,argv);
  parseMMCtrlFile(inputFile);
}

// ----------------------------------------------------------------------
void setCurrentParameters(MoveParameters * params)
{
  // ------------------------------------------------
  // Fill in 'params' with the parsed values
  // ------------------------------------------------
  cout << "Fill in the parameters from the control file...\n";

  // mesh and geometry parameters
  // -----------------------------

  ClassParameter& input = Parser::instance().get("Input"); 
  ClassParameter& inputMesh = input.get("Mesh");
  params->meshFileName = inputMesh.getString("MshFile");
  params->geoFileName  = inputMesh.getString("GeoFile");
  
  // global output parameters
  // ------------------------

  ClassParameter& outputParam = Parser::instance().get("Output");
  params->outType = outputParam.getToken("Type");
  params->outFrequency = outputParam.getInteger("Frequency");
  params->outPrefix = outputParam.getString("Prefix");

  // Global task parameters
  // ----------------------

  ClassParameter& taskParams = Parser::instance().get("Task");
  params->task = taskParams.getToken("Type");

  // debugging parmeters
  // --------------------

  ClassParameter& debugParams = Parser::instance().get("Debug");
  params->debugLevel = debugParams.getInteger("DebugLevel");

  params->openJournal = false;
  string openJournal = debugParams.getToken("OpenJournal");
  if ( !strcmp(openJournal.c_str(),"Yes") ) params->openJournal = true;

  params->referenceJournal = debugParams.getString("ReferenceJournal");

  params->sliverReports = false;
  string sliverReports = debugParams.getToken("ReportSlivers");
  if ( !strcmp(sliverReports.c_str(),"Yes") ) params->sliverReports = true;

  params->testSliverOperators = false;
  string sliverOps = debugParams.getToken("TestSliverOperators");
  if ( !strcmp(sliverOps.c_str(),"Yes") ) params->testSliverOperators = true;

  // Mesh adaptation
  // ----------------

  ClassParameter& adParam = Parser::instance().get("Adaptation");

  params->adapt = true;
  string enable = adParam.getToken("Enable");
  if ( !strcmp(enable.c_str(),"No") ) params->adapt = false;
  params->adaptAlways = true;
  string always = adParam.getToken("AdaptEveryIteration");
  if ( !strcmp(always.c_str(),"No") ) params->adaptAlways = false;
  params->adaptQualityThr = adParam.getDouble("QualityThreshold");

  params->maxInnerIter = adParam.getInteger("MaxNumberOfInnerIteration");
  params->lowerLength = adParam.getDouble("LowerLengthBound");
  params->upperLength = adParam.getDouble("UpperLengthBound");
  params->infLength = adParam.getDouble("InfiniteEdgeLength");
  params->swapMinImproveRatio = adParam.getDouble("SwapMinImproveRatio");
  params->sliverQuality = adParam.getDouble("SliverQuality");

  params->splitCanMakeSlivers = false;
  string slivSpl = adParam.getToken("SplitCanMakeSlivers");
  if ( !strcmp(slivSpl.c_str(),"Yes") ) params->splitCanMakeSlivers = true;
  params->makeSliverInSplitLenSqBound = adParam.getDouble("SliverUpperLengthSqBound");

  params->collapseCanMakeSlivers = false;
  string slivColl = adParam.getToken("CollapseCanMakeSlivers");
  if ( !strcmp(slivColl.c_str(),"Yes") ) params->collapseCanMakeSlivers = true;
  params->makeSliverInCollapseLenSqBound = adParam.getDouble("SliverLowerLengthSqBound");

  params->collapseOnBoundary = false;
  string clpOnB = adParam.getToken("CollapseOnBoundary");
  if ( !strcmp(clpOnB.c_str(),"Yes") ) params->collapseOnBoundary = true;
  params->clpOnBdryTolerance = adParam.get("BoundaryCollapses").getDouble("Tolerance");

  params->swapOnBoundary = false;
  string swapOnB = adParam.getToken("SwapOnBoundary");
  if ( !strcmp(swapOnB.c_str(),"Yes") ) params->swapOnBoundary = true;
  params->swapOnBdryTolerance = adParam.get("BoundarySwaps").getDouble("Tolerance");

  params->trackGeometry = false;
  string trackGeo = adParam.getToken("TrackGeometry");
  if ( !strcmp(trackGeo.c_str(),"Yes") ) {
    if ( params->geoFileName.empty() ) {
      printf("Error: you have to provide a geometry file if you want to track geometry\n");
      throw;
    }
    params->trackGeometry = true;
  }
  ClassParameter& snapParams = adParam.get("VertexSnapping");
  params->snap_force = false;
  if ( !strcmp( snapParams.getToken("ForceRelocation"), "Yes" ) ) {
    params->snap_force = true;
  }
  params->snap_check = false;
  if ( !strcmp( snapParams.getToken("StrictChecking"), "Yes" ) ) {
    params->snap_check = true;
  }
  params->snap_cavityIsMesh = false;
  if ( !strcmp( snapParams.getToken("IsMeshTheCavity"), "Yes" ) ) {
    params->snap_cavityIsMesh = true;
  }
  params->snap_thickness = snapParams.getInteger("CavityThickness");
  params->snap_chi       = snapParams.getDouble("StiffnessAlteration");
  
  params->SFSmoothing = false;
  string smoothSF = adParam.getToken("SmoothSizeField");
  if ( !strcmp(smoothSF.c_str(),"Yes") ) params->SFSmoothing = true;
  ClassParameter& smoothParams = adParam.get("SizeFieldSmoothing");
  params->SFSmoothGrad = smoothParams.getDouble("MaximumGradient");
  
  // Node motion
  // -----------

  ClassParameter& reposParams = Parser::instance().get("NodeMotion");
  params->nodeMotionType = reposParams.getToken("Type");
  ClassParameter& elParams = reposParams.get("ElasticParameters");
  params->elasticStiffnessAlteration = elParams.getDouble("StiffnessAlteration");
  params->elasticSubAdapt = false;
  if ( !strcmp( elParams.getToken("AllowSubAdaptation"), "Yes" ) ) {
    params->elasticSubAdapt = true;
  }
  params->elasticSubAdaptQualThr = elParams.getDouble("SubAdaptationQualityThreshold");
  params->elasticIsMeshTheCavity = false;
  if ( !strcmp( elParams.getToken("IsMeshTheCavity"), "Yes" ) ) {
    params->elasticIsMeshTheCavity = true;
  }
  params->elasticCavityThickness = elParams.getInteger("CavityThickness");

  // Time parameters
  // ---------------

  ClassParameter& temporalParameters = Parser::instance().get("TimeSpecifications");
  params->timeStep       = temporalParameters.getDouble("TimeStep");
  params->finalTime      = temporalParameters.getDouble("FinalTime");
  params->maxNbTimeSteps = temporalParameters.getInteger("MaxTimeSteps");
  
  // Moving objects
  // --------------

  ClassParameter& objsParam = Parser::instance().get("MovingObjects");
  for (int iObj = 0; iObj < objsParam.get("Object").size(); iObj++)
    {
      ClassParameter& objParam = objsParam.get("Object",iObj);

      ObjectDef object;

      // --- Name ---

      object.name = objParam.getTag();

      // --- Geometry ---

      ClassParameter& geoPar = objParam.get("Geometry");
      string geomType = geoPar.getToken("GeometricalType");
      if     (geomType == "Vertex")  object.geomLevel=0;
      else if(geomType == "Edge")    object.geomLevel=1;
      else if(geomType == "Face")    object.geomLevel=2;
      else if(geomType == "Region")  object.geomLevel=3;
      else
        printf("Error: unknown geometrical type %s",geomType.c_str());
      // get tags of entities
      //   vector<int> geomTags = geoPar.getIntegerList("GeometricalTags");
      vector<pair<int,int> > geomTags = geoPar.getIntegerRange("GeometricalTags");
      for (unsigned int iTag = 0; iTag < geomTags.size(); iTag++ )
        for (int iTag2 = (geomTags[iTag]).first; iTag2 <= (geomTags[iTag]).second; iTag2++)
          object.geomTags.insert(iTag2);

      // --- Kinematics ---

      ClassParameter& kinPar = objParam.get("Kinematics");
      string kinType = kinPar.getToken("KinematicsSpecification");
      if ( !strcmp(kinType.c_str(),"Displacement") ) {
        object.kinType = KT_DISPLACEMENT;
        cout << "Displacement formulation for object is wrong "
             << "when using geometry tracking\n";
        ClassParameter& kinParDx = kinPar.get("Displacement");
        object.kinExpression = kinParDx.getStringList("Formulation");
      }

      else if ( !strcmp(kinType.c_str(),"Velocity") ) {
        object.kinType = KT_VELOCITY;
        ClassParameter& kinParV = kinPar.get("Velocity");
        object.kinExpression = kinParV.getStringList("Formulation");
      }
      else throw;

      if (object.kinExpression.size() != 3) {
        cerr << "Error: exactly 3 components are required in the displacement/velocity strings list\n";
        throw;
      }

      // --- Local size fields ---

      for (int iSF = 0; iSF < objParam.get("LocalSizeField").size(); iSF++) {
        ClassParameter& sizeP = objParam.get("LocalSizeField",iSF);
        
        LocalSFDef localSF;

        localSF.name = sizeP.getTag();

        string distComp = sizeP.getToken("DistanceComputation");
        localSF.distToFaces = false;
        if ( !strcmp( distComp.c_str(), "ToWallFaces" ) ) {
          localSF.distToFaces = true;
        }

        localSF.radius = sizeP.getDouble("Radius");

        string iso = sizeP.getToken("Orientation");
        if ( !strcmp(iso.c_str(),"Isotropic") ) {
          localSF.isotropic = true;
          ClassParameter& isoP = sizeP.get("Isotropic");
          localSF.sizeN = isoP.getString("Size");
        }
        else {
          localSF.isotropic = false;
          ClassParameter& anisoP = sizeP.get("Anisotropic");
          localSF.sizeN  = anisoP.getString("NormalSize");
          localSF.sizeT  = anisoP.getString("TangentSize");
        }

        ClassParameter& limiter = sizeP.get("CurvatureBasedLimiter");
        localSF.limit = false;
        if ( !strcmp(limiter.getToken("Enable"), "Yes") ) {
          localSF.limit = true;
          localSF.limitSizeTg = limiter.getDouble("TangentSizeOverRadius");
          localSF.maxCurv     = limiter.getDouble("CurvatureLimiter");
        }
        
        object.sizes.push_back(localSF);
      }

      params->objects.push_back(object);
    }

  // Size Field parameters
  // ---------------------

  for (int iSF = 0; iSF < Parser::instance().get("SizeField").size(); iSF++)
    {
      ClassParameter& sfParam = Parser::instance().get("SizeField",iSF);
      string sfType = sfParam.getToken("Type");

      SizeFieldDef sfDef;

      // Analytical:
      if ( !strcmp(sfType.c_str(),"Analytical") )
        {
          sfDef.type = MAd::ANALYTICALSFIELD;
          
          ClassParameter& anParam = sfParam.get("Analytical");
          string orient = anParam.getToken("OrientationType");
          if ( !strcmp(orient.c_str(),"Isotropic") ) {
            sfDef.orientation = ORT_ISOTROPIC;
            sfDef.isoSize = anParam.get("Isotropic").getString("Length");
          }
          else {
            sfDef.orientation = ORT_ANISOTROPIC;
            ClassParameter& an_anisoParam = anParam.get("Anisotropic");
            
            sfDef.anisoSize = an_anisoParam.getStringList("Lengths");
            sfDef.anisoDir0 = an_anisoParam.getStringList("Direction1");
            sfDef.anisoDir1 = an_anisoParam.getStringList("Direction2");
            sfDef.anisoDir2 = an_anisoParam.getStringList("Direction3");
            
            if ( sfDef.anisoSize.size() != 3 || 
                 sfDef.anisoDir0.size() != 3 || 
                 sfDef.anisoDir1.size() != 3 ||  
                 sfDef.anisoDir2.size() != 3 ) {
              cerr << "Invalid number of strings in one of the lists (anisotropic case)\n";
              throw;
            }
          }
        }

      // Piecewise linear:
      else if ( !strcmp(sfType.c_str(),"PiecewiseLinear") )
        {
          sfDef.type = MAd::DISCRETESFIELD;
          ClassParameter& sf_pwl = sfParam.get("PiecewiseLinear");
          sfDef.pwlSource = sf_pwl.getToken("Source");
          ClassParameter& sf_pwl_curv = sf_pwl.get("Curvature");
  
          sfDef.curv_aniso = false;
          if ( !strcmp( sf_pwl_curv.getToken("Anisotropic"), "Yes") ) {
            sfDef.curv_aniso = true;
          }
          sfDef.curv_alpha = sf_pwl_curv.getDouble("Alpha");
          sfDef.curv_hMin  = sf_pwl_curv.getDouble("MinimalLength");
        }

      // background:
      else if ( !strcmp(sfType.c_str(),"Background") )
        {
          sfDef.type = MAd::BACKGROUNDSFIELD;

          ClassParameter& sf_bg = sfParam.get("Background");
          sfDef.mshFile = sf_bg.getString("FileName");
          
        }

      else throw;

      params->sizes.push_back(sfDef);
    }
}

// -----------------------------------------------------------------------------

#endif

#endif
