#include <vtkUnstructuredGrid.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <set>
#include <string>
#include "MeshGeneration/meshGen.H"
#include "MeshGeneration/snappymeshGen.H"
#include "MeshGeneration/snappymeshParams.H"

// openfoam headers
#include <IOmanip.H>
#include <MeshedSurface.H>
#include <Time.H>
#include <UnsortedMeshedSurface.H>
#include <argList.H>
#include <cellModeller.H>
#include <decompositionMethod.H>
#include <faceSet.H>
#include <fileName.H>
#include <fvCFD.H>
#include <fvMesh.H>
#include <fvMeshDistribute.H>
#include <fvMeshTools.H>
#include <globalIndex.H>
#include <meshRefinement.H>
#include <motionSmoother.H>
#include <noDecomp.H>
#include <polyTopoChange.H>
#include <refinementParameters.H>
#include <surfZoneIdentifierList.H>
#include <uindirectPrimitivePatch.H>
#include <vtkSetWriter.H>
#include <wallPolyPatch.H>
#include <foamVTKTopo.H>
#include <decompositionModel.H>
#include <processorMeshes.H>
#include <profiling.H>
#include <snappyVoxelMeshDriver.H>


// snappyHexMesh Headers
#include <layerParameters.H>
#include <refinementFeatures.H>
#include <refinementSurfaces.H>
#include <searchableSurfaces.H>
#include <shellSurfaces.H>
#include <snapParameters.H>
#include <snappyLayerDriver.H>
#include <snappyRefineDriver.H>
#include <snappySnapDriver.H>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// This class has all the methods from snappyHexMesh. Currently, Surface
// Simplyfy method is not enabled.

snappymeshGen::snappymeshGen()  // default constructor definition
{
  // Default Meshing Parameters
  _params = new snappymeshParams();
  defaults = true;

  // Initialization Tasks
  initialize();
}

// constructor definition with parameters input
snappymeshGen::snappymeshGen(snappymeshParams *params)
    : defaults(false), _params(params) {
  initialize();  // Foam Initialization
}

snappymeshGen::~snappymeshGen()  // destructor definition
{
  if (defaults) delete _params;
}

void snappymeshGen::initialize() {
  // create dictionaries needed
  createControlDict();
  createSnappyDict();
  createfvSchemesDict();
  createfvSolutionDict();
}

void snappymeshGen::createControlDict() {
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating system directory for the cfMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path) + "/controlDict");
  std::string contText =
      "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    location  \"system\";\n\
    object    controlDict;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
deltaT  1;\n\n\
startTime   0;\n\n\
writeInterval   1;\n\n\
// ***************************************************************** //";
  contDict << contText;
  contDict.close();
}

void snappymeshGen::createfvSchemesDict() {
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating system directory for the cfMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path) + "/fvSchemes");
  std::string contText =
      "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    object    fvSchemes;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
gradSchemes\n\
{\n\
    default         Gauss linear;\n\
    grad(p)         Gauss linear;\n\
}\n\
\n\
divSchemes\n\
{\n\
    default         none;\n\
    div(phi,U)      Gauss linear;\n\
}\n\
\n\
laplacianSchemes\n\
{\n\
    default         none;\n\
    laplacian(nu,U) Gauss linear corrected;\n\
    laplacian((1|A(U)),p) Gauss linear corrected;\n\
}\n\
// ******************************************************************* //";
  contDict << contText;
  contDict.close();
}

void snappymeshGen::createfvSolutionDict() {
  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating system directory for the cfMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  std::ofstream contDict;
  contDict.open(std::string(dir_path) + "/fvSolution");
  std::string contText =
      "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: cfMesh interface                      |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    object    fvSolution;\n\
}\n\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
\n\
// ********************************************************************** //";
  contDict << contText;
  contDict.close();
}

void snappymeshGen::createSnappyDict() {
  // creating snappyHexMeshDict

  // creating a base system directory
  const char dir_path[] = "./system";
  boost::filesystem::path dir(dir_path);
  try {
    boost::filesystem::create_directory(dir);
  } catch (boost::filesystem::filesystem_error &e) {
    std::cerr << "Problem in creating system directory for the snappyHexMesh"
              << "\n";
    std::cerr << e.what() << std::endl;
    throw;
  }

  // creating mesh dictionary file
  std::ofstream contDict;
  contDict.open(std::string(dir_path) + "/snappyHexMeshDict");
  // header
  std::string contText =
      "\
/*--------------------------------*- C++ -*----------------------------------*\n\
| =========                 |                                                |\n\
| \\\\      /  F ield         | NEMoSys: snappyHexMesh interface               |\n\
|  \\\\    /   O peration     |                                                |\n\
|   \\\\  /    A nd           |                                                |\n\
|    \\\\/     M anipulation  |                                                |\n\
\\*---------------------------------------------------------------------------*/\n\
\n\
FoamFile\n\
{\n\
    version   2.0;\n\
    format    ascii;\n\
    class     dictionary;\n\
    location  \"system\";\n\
    object    snappyHexMeshDict;\n\
}\n\n";

  contText =
      contText +
      "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";

  // initial booleans
  if ((_params->withCastMesh)) {
    contText = contText + "\ncastellatedMesh  true;\n";
  } else {
    contText = contText + "\ncastellatedMesh  false;\n";
  }

  if ((_params->withSnap)) {
    contText = contText + "\nsnap       true;\n";
  } else {
    contText = contText + "\nsnap       false;\n";
  }

  if ((_params->withLayers)) {
    contText = contText + "\naddLayers      true;\n";
  } else {
    contText = contText + "\naddLayers      false;\n";
  }

  // Geometry declaration
  contText = contText + "\ngeometry\n{\n";
  contText = contText + "\t" + (_params->geomFileName);
  contText = contText + "\n\t{\n\t\ttype\ttriSurfaceMesh;\n";

  if (_params->geoDef.withMultiPatches) {
    if ((_params->geoDef.stlPatchDefs).size() == 0) {
      std::cerr << "Error reading in custom patches from JSON input!"
                << std::endl;
      throw;
    }
    contText = contText + "\t\tregions\n\t\t{\n";
    for (auto pt = (_params->geoDef.stlPatchDefs).begin();
         pt != (_params->geoDef.stlPatchDefs).end(); pt++) {
      contText = contText + "\t\t\t" + pt->STLPatchName + "  \t{name " +
                 pt->snappyPatchName + ";}\n";
    }
    contText = contText + "\t\t}\n";
  } else {
    contText = contText + "\t\tname\t" + (_params->geoDef.singleSolidPatch) + ";\n";
  }

  // Searchable shapes go here
  if ((_params->geoDef.srchShape).size() > 0) {
    for (const auto &shape : _params->geoDef.srchShape) {
      contText = contText + "\t\t" + (shape->patchName);
      contText = contText + "\n\t\t{";
      contText = contText + "\n\t\t\ttype " + (shape->getType()) + ";\n";

      if (auto box = std::dynamic_pointer_cast<shmSearchableBox>(shape)) {
        contText = contText + "\t\t\tmin (" +
                   std::to_string(box->minBound.at(0)) + " " +
                   std::to_string(box->minBound.at(1)) + " " +
                   std::to_string(box->minBound.at(2)) + ");\n";
        contText = contText + "\t\t\tmax (" +
                   std::to_string(box->maxBound.at(0)) + " " +
                   std::to_string(box->maxBound.at(1)) + " " +
                   std::to_string(box->maxBound.at(2)) + ");\n";
      } else if (auto sphere =
                     std::dynamic_pointer_cast<shmSearchableSphere>(shape)) {
        contText = contText + "\t\t\tcentre (" +
                   std::to_string(sphere->center.at(0)) + " " +
                   std::to_string(sphere->center.at(1)) + " " +
                   std::to_string(sphere->center.at(2)) + ");\n";
        contText = contText + "\t\t\tradius " +
                   std::to_string(sphere->radius) + ";\n";
      } else if (auto cyl =
                     std::dynamic_pointer_cast<shmSearchableCylinder>(shape)) {
        contText = contText + "\t\t\tpoint1 (" +
                   std::to_string(cyl->axisPoint1.at(0)) + " " +
                   std::to_string(cyl->axisPoint1.at(1)) + " " +
                   std::to_string(cyl->axisPoint1.at(2)) + ");\n";
        contText = contText + "\t\t\tpoint2 (" +
                   std::to_string(cyl->axisPoint2.at(0)) + " " +
                   std::to_string(cyl->axisPoint2.at(1)) + " " +
                   std::to_string(cyl->axisPoint2.at(2)) + ");\n";
        contText =
            contText + "\t\t\tradius " + std::to_string(cyl->radius) + ";\n";
      }
      contText = contText + "\n\t\t}\n";
    }
  }
  contText = contText + "\t}\n";
  contText = contText + "};\n";
  // End of geometry definition

  // Castellated Mesh Controls
  contText = contText + "\n\ncastellatedMeshControls\n{\n";
  contText = contText + "\tmaxLocalCells\t" +
             std::to_string(_params->castMeshControls.maxLCells) + ";\n";
  contText = contText + "\tmaxGlobalCells\t" +
             std::to_string(_params->castMeshControls.maxGCells) + ";\n";
  contText = contText + "\tnCellsBetweenLevels\t" +
             std::to_string(_params->castMeshControls.cellsBetnLvls) + ";\n";
  contText = contText + "\tminRefinementCells\t" +
             std::to_string(_params->castMeshControls.minRefCells) + ";\n";

  // Features using emesh file
  contText = contText + "\n\tfeatures\n\t(\n";
  if ((_params->castMeshControls.ftrEdge).size() > 0) {
    for (auto pt = (_params->castMeshControls.ftrEdge).begin();
         pt != (_params->castMeshControls.ftrEdge).end(); pt++) {
      contText = contText + "\t\t{\n";
      contText = contText + "\t\t\tfile \"" + pt->fileName + "\";\n";
      contText = contText + "\t\t\tlevels ((" + std::to_string(pt->minLvl) +
                 " " + std::to_string(pt->maxLvl) + "));\n";
      contText = contText + "\t\t}\n";
    }
  }
  contText = contText + "\t);\n";

  contText = contText + "\n\trefinementSurfaces\n\t{\n";
  if (_params->geoDef.withMultiPatches)
    contText = contText + "\t\t" + (_params->geomFileName) + "\n\t\t{";
  else
    contText =
        contText + "\t\t" + (_params->geoDef.singleSolidPatch) + "\n\t\t{";
  contText = contText + "\n\t\t\tlevel (" +
             std::to_string(_params->castMeshControls.refSurfLvlMin) + " " +
             std::to_string(_params->castMeshControls.refSurfLvlMax) + ");\n";

  if ((_params->castMeshControls.withCellZones)) {
    contText = contText + "\t\t\tfaceZone\t" +
               (_params->geoDef.singleSolidPatch) + ";\n";
    contText = contText + "\t\t\tcellZone\t" +
               (_params->geoDef.singleSolidPatch) + ";\n";
    contText = contText + "\t\t\tcellZoneInside\tinside;\n";
  }

  if ((_params->geoDef.withMultiPatches)) {
    if (((_params->castMeshControls.surfRefs).size() > 0)) {
      contText = contText + "\n\t\t\tregions\n";
      contText = contText + "\t\t\t{";

      for (auto pt = (_params->castMeshControls.surfRefs).begin();
           pt != (_params->castMeshControls.surfRefs).end(); pt++) {
        if (pt->patchType == "NO") {
          contText = contText + "\n\t\t\t\t" + (pt->refPatchName) +
                     "\t{level (" + std::to_string((int)pt->minLvl) + " " +
                     std::to_string(pt->maxLvl) +
                     "); patchInfo { type patch; }}";
        } else {
          contText = contText + "\n\t\t\t\t" + (pt->refPatchName) +
                     "\t{level (" + std::to_string((int)pt->minLvl) + " " +
                     std::to_string(pt->maxLvl) + "); patchInfo { type " +
                     pt->patchType + "; }}";
        }
      }

      contText = contText + "\n\t\t\t}\n";
    }
  }

  contText = contText + "\t\t}\n";

  contText = contText + "\t\tgapLevelIncrement " +
             std::to_string(_params->castMeshControls.castMeshGpLvl) + +";\n";

  contText = contText + "\t}\n";

  contText = contText + "\tresolveFeatureAngle\t" +
             std::to_string(_params->castMeshControls.featAngle) + ";\n";

  contText = contText + "\tgapLevelIncrement\t" +
             std::to_string(_params->castMeshControls.gapLvlInc) + ";\n";

  contText = contText + "\n\tplanarAngle " +
             std::to_string(_params->castMeshControls.planarAngle) + ";\n";

  contText = contText + "\n\trefinementRegions\n\t{\n";

  if (((_params->castMeshControls.geomRefs).size() > 0)) {
    for (auto pt = (_params->castMeshControls.geomRefs).begin();
         pt != (_params->castMeshControls.geomRefs).end(); pt++) {
      contText = contText + "\n\t\t" + (pt->patchName);
      contText = contText + "\n\t\t{";
      contText = contText + "\n\t\t\tmode " + (pt->mode) + ";\n";
      contText = contText + "\t\t\tlevels ((" + std::to_string(pt->minLvl) +
                 " " + std::to_string(pt->maxLvl) + "));\n";
      contText = contText + "\t\t}\n";
    }
  }

  contText = contText + "\t}";

  contText = contText + "\n\tlocationInMesh\t(" +
             std::to_string(_params->castMeshControls.locMesh.at(0)) + " " +
             std::to_string(_params->castMeshControls.locMesh.at(1)) + " " +
             std::to_string(_params->castMeshControls.locMesh.at(2)) + ");\n";
  if ((_params->castMeshControls.alwFreeZone) == 1)
    contText = contText + "\tallowFreeStandingZoneFaces\ttrue;\n";
  if ((_params->castMeshControls.alwFreeZone) == 0)
    contText = contText + "\tallowFreeStandingZoneFaces\tfalse;\n";
  contText = contText + "}\n";

  // Snap Controls
  contText = contText + "\n\nsnapControls\n{\n";
  contText = contText + "\tnSmoothPatch\t" +
             std::to_string(_params->snapControls.nSmoothPatch) + ";\n";
  contText = contText + "\ttolerance\t" +
             std::to_string(_params->snapControls.tolerance) + ";\n";
  contText = contText + "\tnSolveIter\t" +
             std::to_string(_params->snapControls.nSolveIter) + ";\n";
  contText = contText + "\tnRelaxIter\t" +
             std::to_string(_params->snapControls.nRelaxIter) + ";\n";

  // New Features
  contText = contText + "\tnFeatureSnapIter\t" +
             std::to_string(_params->snapControls.nFeatureSnapIter) + ";\n";

  if (_params->snapControls.implicitFeatureSnap)
    contText = contText + "\timplicitFeatureSnap\t true;\n";
  else
    contText = contText + "\timplicitFeatureSnap\t false;\n";

  if (_params->snapControls.explicitFeatureSnap)
    contText = contText + "\texplicitFeatureSnap\t true;\n";
  else
    contText = contText + "\texplicitFeatureSnap\t false;\n";

  if (_params->snapControls.multiRegionFeatureSnap)
    contText = contText + "\tmultiRegionFeatureSnap\t true;\n";
  else
    contText = contText + "\tmultiRegionFeatureSnap\t false;\n";

  contText = contText + "}\n";

  // Layer Controls
  contText = contText + "\n\naddLayersControls\n{\n";
  if ((_params->layerControls.relSize) == 1)
    contText = contText + "\trelativeSizes\ttrue;\n";
  if ((_params->layerControls.relSize) == 0)
    contText = contText + "\trelativeSizes\tfalse;\n";

  // Layers
  contText = contText + "\tlayers\n\t{\n";
  if ((_params->layerControls.layerVec).size() > 0) {
    for (auto pt = (_params->layerControls.layerVec).begin();
         pt != (_params->layerControls.layerVec).end(); pt++) {
      contText = contText + "\t\t" + pt->patchName + "\n\t\t{\n";
      contText = contText + "\t\t\tnSurfaceLayers " +
                 std::to_string(pt->nSurfaceLayers) + ";\n";
      contText = contText + "\t\t\texpansionRatio " +
                 std::to_string(pt->expansionRatio) + ";\n";
      contText = contText + "\t\t\tfinalLayerThickness " +
                 std::to_string(pt->finalLayerThickness) + ";\n";
      if (pt->firstLyrThickness > 0)
        contText = contText + "\t\t\tfirstLayerThickness " +
                   std::to_string(pt->firstLyrThickness) + ";\n";
      if (pt->thickness > 0)
        contText = contText + "\t\t\tthickness " +
                   std::to_string(pt->thickness) + ";\n";
      contText = contText + "\t\t\tminThickness " +
                 std::to_string(pt->minThickness) + ";\n";
      contText = contText + "\t\t}\n";
    }
  }
  contText = contText + "\t}\n";

  contText = contText + "\texpansionRatio\t" +
             std::to_string(_params->layerControls.expRatio) + ";\n";
  contText = contText + "\tfinalLayerThickness\t" +
             std::to_string(_params->layerControls.finLThick) + ";\n";
  contText = contText + "\tminThickness\t" +
             std::to_string(_params->layerControls.minThick) + ";\n";
  if (_params->layerControls.firstLyrThickness > 0)
    contText = contText + "\tfirstLayerThickness\t" +
               std::to_string(_params->layerControls.firstLyrThickness) + ";\n";
  if (_params->layerControls.thickness > 0)
    contText = contText + "\tthickness\t" +
               std::to_string(_params->layerControls.thickness) + ";\n";
  contText = contText + "\tnGrow\t" +
             std::to_string(_params->layerControls.nGrow) + ";\n";
  contText = contText + "\tfeatureAngle\t" +
             std::to_string(_params->layerControls.featAngle) + ";\n";
  contText = contText + "\tslipFeatureAngle\t" +
             std::to_string(_params->layerControls.slipFeatureAngle) + ";\n";
  contText = contText + "\tnRelaxIter\t" +
             std::to_string(_params->layerControls.relaxIter) + ";\n";
  contText = contText + "\tnSmoothSurfaceNormals\t" +
             std::to_string(_params->layerControls.smthSurfNorm) + ";\n";
  contText = contText + "\tnSmoothNormals\t" +
             std::to_string(_params->layerControls.smthNorm) + ";\n";
  contText = contText + "\tnSmoothThickness\t" +
             std::to_string(_params->layerControls.smthThick) + ";\n";
  contText = contText + "\tmaxFaceThicknessRatio\t" +
             std::to_string(_params->layerControls.maxFcTR) + ";\n";
  contText = contText + "\tmaxThicknessToMedialRatio\t" +
             std::to_string(_params->layerControls.maxThickTMR) + ";\n";
  contText = contText + "\tminMedialAxisAngle\t" +
             std::to_string(_params->layerControls.minMedAngl) + ";\n";
  contText = contText + "\tnBufferCellsNoExtrude\t" +
             std::to_string(_params->layerControls.bufferCells) + ";\n";
  contText = contText + "\tnLayerIter\t" +
             std::to_string(_params->layerControls.nIter) + ";\n";
  contText = contText + "\tnRelaxedIter\t" +
             std::to_string(_params->layerControls.nRelaxedIter) + ";\n";
  if (_params->layerControls.nMedialAxisIter > 0)
    contText = contText + "\tnMedialAxisIter\t" +
               std::to_string(_params->layerControls.nMedialAxisIter) + ";\n";
  if (_params->layerControls.nSmoothDisplacement > 0)
    contText = contText + "\tnSmoothDisplacement\t" +
               std::to_string(_params->layerControls.nSmoothDisplacement) +
               ";\n";
  contText = contText + "}\n";

  // Mesh Quality Controls
  std::ostringstream minimumVol;
  minimumVol << _params->qualityControls.minVol;
  std::ostringstream minimumTetQuality;
  minimumTetQuality << _params->qualityControls.minTetQ;

  contText = contText + "\n\nmeshQualityControls\n{\n";
  contText = contText + "\tmaxNonOrtho\t" +
             std::to_string(_params->qualityControls.maxNonOrtho) + ";\n";
  contText = contText + "\tmaxBoundarySkewness\t" +
             std::to_string(_params->qualityControls.maxBndrySkew) + ";\n";
  contText = contText + "\tmaxInternalSkewness\t" +
             std::to_string(_params->qualityControls.maxIntSkew) + ";\n";
  contText = contText + "\tmaxConcave\t" +
             std::to_string(_params->qualityControls.maxConc) + ";\n";
  contText = contText + "\tminVol\t" + minimumVol.str() + ";\n";
  contText = contText + "\tminTetQuality\t" + minimumTetQuality.str() + ";\n";
  contText = contText + "\tminArea\t" +
             std::to_string(_params->qualityControls.minArea) + ";\n";
  contText = contText + "\tminTwist\t" +
             std::to_string(_params->qualityControls.minTwist) + ";\n";
  contText = contText + "\tminFaceWeight\t" +
             std::to_string(_params->qualityControls.minFaceW) + ";\n";
  contText = contText + "\tminVolRatio\t" +
             std::to_string(_params->qualityControls.minVolRto) + ";\n";
  contText = contText + "\tminDeterminant\t" +
             std::to_string(_params->qualityControls.minDet) + ";\n";
  contText = contText + "\tminTriangleTwist\t" +
             std::to_string(_params->qualityControls.minTriTwist) + ";\n";
  contText = contText + "\tnSmoothScale\t" +
             std::to_string(_params->qualityControls.smoothScale) + ";\n";
  contText = contText + "\terrorReduction\t" +
             std::to_string(_params->qualityControls.errReduction) + ";\n";
  contText = contText + "\trelaxed\n\t{\n";
  contText = contText + "\t\tmaxNonOrtho \t " +
             std::to_string(_params->qualityControls.maxNonOrtho) + ";\n\t}\n";
  contText = contText + "}\n";

  contText = contText + "\nmergeTolerance\t" +
             std::to_string(_params->mergeTol) + ";\n\n";

  contText = contText + "writeFlags\n(\n";
  contText = contText + "\tscalarLevels\n";
  contText = contText + "\tlayerSets\n";
  contText = contText + "\tlayerFields\n);\n";

  contText =
      contText +
      "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
  contDict << contText;
  contDict.close();
}

// Convert size (as fraction of defaultCellSize) to refinement level
Foam::label snappymeshGen::sizeCoeffToRefinement(
    const Foam::scalar level0Coeff,  // ratio of hex cell size
                                     // v.s. defaultCellSize
    const Foam::scalar sizeCoeff) {
  using namespace Foam;
  return round(::log(level0Coeff / sizeCoeff) / ::log(2));
}

Foam::autoPtr<refinementSurfaces> snappymeshGen::createRefinementSurfaces(
    const Foam::searchableSurfaces &allGeometry,
    const Foam::dictionary &surfacesDict,
    const Foam::dictionary &shapeControlDict,
    const Foam::label gapLevelIncrement, const Foam::scalar level0Coeff) {
  using namespace Foam;

  autoPtr<refinementSurfaces> surfacePtr;

  // Count number of surfaces.
  label surfI = 0;
  forAll (allGeometry.names(), geomI) {
    const word &geomName = allGeometry.names()[geomI];

    if (surfacesDict.found(geomName)) {
      surfI++;
    }
  }

  labelList surfaces(surfI);
  wordList names(surfI);
  PtrList<surfaceZonesInfo> surfZones(surfI);

  labelList regionOffset(surfI);

  labelList globalMinLevel(surfI, Zero);
  labelList globalMaxLevel(surfI, Zero);
  labelList globalLevelIncr(surfI, Zero);
  PtrList<dictionary> globalPatchInfo(surfI);
  List<Map<label>> regionMinLevel(surfI);
  List<Map<label>> regionMaxLevel(surfI);
  List<Map<label>> regionLevelIncr(surfI);
  List<Map<scalar>> regionAngle(surfI);
  List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);

  HashSet<word> unmatchedKeys(surfacesDict.toc());

  surfI = 0;
  forAll (allGeometry.names(), geomI) {
    const word &geomName = allGeometry.names()[geomI];

    const entry *ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);

    if (ePtr) {
      const dictionary &shapeDict = ePtr->dict();
      unmatchedKeys.erase(ePtr->keyword());

      names[surfI] = geomName;
      surfaces[surfI] = geomI;

      const searchableSurface &surface = allGeometry[geomI];

      // Find the index in shapeControlDict
      // Invert surfaceCellSize to get the refinementLevel
      const word scsFuncName =
                shapeDict.get<word>("surfaceCellSizeFunction");
      const dictionary& scsDict =
                shapeDict.optionalSubDict(scsFuncName + "Coeffs");
      const scalar surfaceCellSize =
          scsDict.get<scalar>("surfaceCellSizeCoeff");

      const label refLevel =
          sizeCoeffToRefinement(level0Coeff, surfaceCellSize);

      globalMinLevel[surfI] = refLevel;
      globalMaxLevel[surfI] = refLevel;
      globalLevelIncr[surfI] = gapLevelIncrement;

      // Surface zones    
      surfZones.set(surfI,
                    new surfaceZonesInfo(surface,shapeDict,
                    allGeometry.regionNames()[surfaces[surfI]]));

      // Global perpendicular angle
      if (shapeDict.found("patchInfo")) {
        globalPatchInfo.set(surfI, shapeDict.subDict("patchInfo").clone());
      }

      // Per region override of patchInfo

      if (shapeDict.found("regions")) {
        const dictionary &regionsDict = shapeDict.subDict("regions");
        const wordList &regionNames = allGeometry[surfaces[surfI]].regions();

        forAll (regionNames, regionI) {
          if (regionsDict.found(regionNames[regionI])) {
            // Get the dictionary for region
            const dictionary &regionDict =
                regionsDict.subDict(regionNames[regionI]);

            if (regionDict.found("patchInfo")) {
              regionPatchInfo[surfI].insert(
                  regionI, regionDict.subDict("patchInfo").clone());
            }
          }
        }
      }

      // Per region override of cellSize
      if (shapeDict.found("regions")) {
        const dictionary &shapeControlRegionsDict =
            shapeDict.subDict("regions");
        const wordList &regionNames = allGeometry[surfaces[surfI]].regions();

        forAll (regionNames, regionI) {
          if (shapeControlRegionsDict.found(regionNames[regionI])) {
            const dictionary &shapeControlRegionDict =
                shapeControlRegionsDict.subDict(regionNames[regionI]);            
            const word scsFuncName =
                shapeDict.get<word>("surfaceCellSizeFunction");
            const dictionary &scsDict =
                shapeControlRegionDict.subDict(scsFuncName + "Coeffs");
            const scalar surfaceCellSize =
                            scsDict.get<scalar>("surfaceCellSizeCoeff");

            const label refLevel =
                sizeCoeffToRefinement(level0Coeff, surfaceCellSize);

            regionMinLevel[surfI].insert(regionI, refLevel);
            regionMaxLevel[surfI].insert(regionI, refLevel);
            regionLevelIncr[surfI].insert(regionI, 0);
          }
        }
      }

      surfI++;
    }
  }

  // Calculate local to global region offset
  label nRegions = 0;

  forAll (surfaces, surfI) {
    regionOffset[surfI] = nRegions;
    nRegions += allGeometry[surfaces[surfI]].regions().size();
  }

  // Rework surface specific information into information per global region
  labelList minLevel(nRegions, Zero);
  labelList maxLevel(nRegions, Zero);
  labelList gapLevel(nRegions, -1);
  PtrList<dictionary> patchInfo(nRegions);

  forAll (globalMinLevel, surfI) {
    label nRegions = allGeometry[surfaces[surfI]].regions().size();

    // Initialise to global (i.e. per surface)
    for (label i = 0; i < nRegions; i++) {
      label globalRegionI = regionOffset[surfI] + i;
      minLevel[globalRegionI] = globalMinLevel[surfI];
      maxLevel[globalRegionI] = globalMaxLevel[surfI];
      gapLevel[globalRegionI] =
          maxLevel[globalRegionI] + globalLevelIncr[surfI];

      if (globalPatchInfo.set(surfI)) {
        patchInfo.set(globalRegionI, globalPatchInfo[surfI].clone());
      }
    }

    // Overwrite with region specific information
    forAllConstIter(Map<label>, regionMinLevel[surfI], iter) {
      label globalRegionI = regionOffset[surfI] + iter.key();

      minLevel[globalRegionI] = iter();
      maxLevel[globalRegionI] = regionMaxLevel[surfI][iter.key()];
      gapLevel[globalRegionI] =
          maxLevel[globalRegionI] + regionLevelIncr[surfI][iter.key()];
    }

    const Map<autoPtr<dictionary>> &localInfo = regionPatchInfo[surfI];
    forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter) {
      label globalRegionI = regionOffset[surfI] + iter.key();
      patchInfo.set(globalRegionI, iter()().clone());
    }
  }

  surfacePtr.reset
  (
    new refinementSurfaces
      (
      allGeometry,
      surfaces,
      names,
      surfZones,
      regionOffset,
      minLevel,
      maxLevel,
      gapLevel,
      scalarField(nRegions, -GREAT),  //perpendicularAngle,
      patchInfo,
      false                           //dryRun
    )
  );

  const refinementSurfaces &rf = surfacePtr();

  // Determine maximum region name length
  label maxLen = 0;
  forAll (rf.surfaces(), surfI) {
    label geomI = rf.surfaces()[surfI];
    const wordList &regionNames = allGeometry.regionNames()[geomI];
    forAll (regionNames, regionI) {
      maxLen = Foam::max(maxLen, label(regionNames[regionI].size()));
    }
  }

  Info << Foam::setw(maxLen) << "Region" << Foam::setw(10) << "Min Level"
       << Foam::setw(10) << "Max Level" << Foam::setw(10) << "Gap Level" << nl
       << Foam::setw(maxLen) << "------" << Foam::setw(10) << "---------"
       << Foam::setw(10) << "---------" << Foam::setw(10) << "---------"
       << endl;

  forAll (rf.surfaces(), surfI) {
    label geomI = rf.surfaces()[surfI];

    Info << rf.names()[surfI] << ':' << nl;

    const wordList &regionNames = allGeometry.regionNames()[geomI];

    forAll (regionNames, regionI) {
      label globalI = rf.globalRegion(surfI, regionI);

      Info << Foam::setw(maxLen) << regionNames[regionI] << Foam::setw(10)
           << rf.minLevel()[globalI] << Foam::setw(10) << rf.maxLevel()[globalI]
           << Foam::setw(10) << rf.gapLevel()[globalI] << endl;
    }
  }

  return surfacePtr;
}

void snappymeshGen::writeMesh(
    const Foam::string &msg, const Foam::meshRefinement &meshRefiner,
    const Foam::meshRefinement::debugType debugLevel,
    const Foam::meshRefinement::writeType writeLevel) {
  using namespace Foam;
  const fvMesh &mesh = meshRefiner.mesh();

  meshRefiner.printMeshInfo(debugLevel, msg);
  Info << "Writing mesh to time " << meshRefiner.timeName() << endl;

  processorMeshes::removeFiles(mesh);
    if (!debugLevel && !(writeLevel&meshRefinement::WRITELAYERSETS))
      topoSet::removeFiles(mesh);
    refinementHistory::removeFiles(mesh);

  meshRefiner.write(
      debugLevel,
      meshRefinement::writeType(writeLevel | meshRefinement::WRITEMESH),
      mesh.time().path() / meshRefiner.timeName());
  Info << "Wrote mesh in = " << mesh.time().cpuTimeIncrement() << " s." << endl;
}

void snappymeshGen::removeZeroSizedPatches(Foam::fvMesh &mesh) {
  using namespace Foam;
  // Remove any zero-sized ones. Assumes
  // - processor patches are already only there if needed
  // - all other patches are available on all processors
  // - but coupled ones might still be needed, even if zero-size
  //   (e.g. processorCyclic)
  // See also logic in createPatch.
  const polyBoundaryMesh &pbm = mesh.boundaryMesh();

  labelList oldToNew(pbm.size(), -1);
  label newPatchi = 0;
  forAll (pbm, patchi) {
    const polyPatch &pp = pbm[patchi];

    if (!isA<processorPolyPatch>(pp)) {
      if (isA<coupledPolyPatch>(pp) ||
          returnReduce(pp.size(), sumOp<label>())) {
        // Coupled (and unknown size) or uncoupled and used
        oldToNew[patchi] = newPatchi++;
      }
    }
  }

  forAll (pbm, patchi) {
    const polyPatch &pp = pbm[patchi];

    if (isA<processorPolyPatch>(pp)) {
      oldToNew[patchi] = newPatchi++;
    }
  }

  const label nKeepPatches = newPatchi;

  // Shuffle unused ones to end
  if (nKeepPatches != pbm.size()) {
    Info << endl << "Removing zero-sized patches:" << endl << incrIndent;

    forAll (oldToNew, patchi) {
      if (oldToNew[patchi] == -1) {
        Info << indent << pbm[patchi].name() << " type " << pbm[patchi].type()
             << " at position " << patchi << endl;
        oldToNew[patchi] = newPatchi++;
      }
    }
    Info << decrIndent;

    fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
    Info << endl;
  }
}

Foam::scalar snappymeshGen::getMergeDistance(const Foam::polyMesh &mesh,
                                             const Foam::scalar mergeTol,
                                             const bool dryRun) {
  using namespace Foam;
  const boundBox &meshBb = mesh.bounds();
  scalar mergeDist = mergeTol * meshBb.mag();

  Info << nl << "Overall mesh bounding box  : " << meshBb << nl
       << "Relative tolerance         : " << mergeTol << nl
       << "Absolute matching distance : " << mergeDist << nl << endl;

  // check writing tolerance
  if (mesh.time().writeFormat() == IOstream::ASCII && !dryRun) {
    const scalar writeTol =
        std::pow(scalar(10.0), -scalar(IOstream::defaultPrecision()));

    if (mergeTol < writeTol) {
      FatalErrorInFunction
          << "Your current settings specify ASCII writing with "
          << IOstream::defaultPrecision() << " digits precision." << nl
          << "Your merging tolerance (" << mergeTol << ") is finer than this."
          << nl << "Change to binary writeFormat, "
          << "or increase the writePrecision" << endl
          << "or adjust the merge tolerance (mergeTol)." << exit(FatalError);
    }
  }

  return mergeDist;
}

void snappymeshGen::extractSurface(const Foam::polyMesh &mesh,
                                   const Foam::Time &runTime,
                                   const Foam::labelHashSet &includePatches,
                                   const Foam::fileName &outFileName) {
  using namespace Foam;
  const polyBoundaryMesh &bMesh = mesh.boundaryMesh();

  // Collect sizes. Hash on names to handle local-only patches (e.g.
  //  processor patches)
  HashTable<label> patchSize(1024);
  label nFaces = 0;
  for (const label patchi : includePatches) {
    const polyPatch& pp = bMesh[patchi];
    patchSize.insert(pp.name(), pp.size());
    nFaces += pp.size();
  }
  Pstream::mapCombineGather(patchSize, plusEqOp<label>());

  // Allocate zone/patch for all patches
  HashTable<label> compactZoneID(1024);
  forAllConstIters(patchSize, iter) {
    label sz = compactZoneID.size();
    compactZoneID.insert(iter.key(), sz);
  }
  Pstream::mapCombineScatter(compactZoneID);

  // Rework HashTable into labelList just for speed of conversion
  labelList patchToCompactZone(bMesh.size(), -1);
  forAllConstIters(compactZoneID, iter) {
    label patchi = bMesh.findPatchID(iter.key());
    if (patchi != -1) {
      patchToCompactZone[patchi] = iter();
    }
  }

  // Collect faces on zones
  DynamicList<label> faceLabels(nFaces);
  DynamicList<label> compactZones(nFaces);
  for (const label patchi : includePatches) {
    const polyPatch& pp = bMesh[patchi];
    forAll(pp, i) {
      faceLabels.append(pp.start()+i);
      compactZones.append(patchToCompactZone[pp.index()]);
    }
  }

  // Addressing engine for all faces
  uindirectPrimitivePatch allBoundary(
      UIndirectList<face>(mesh.faces(), faceLabels), mesh.points());

  // Find correspondence to master points
  labelList pointToGlobal;
  labelList uniqueMeshPoints;
  autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints(
      allBoundary.meshPoints(), allBoundary.meshPointMap(), pointToGlobal,
      uniqueMeshPoints);

  // Gather all unique points on master
  List<pointField> gatheredPoints(Pstream::nProcs());
  gatheredPoints[Pstream::myProcNo()] =
      pointField(mesh.points(), uniqueMeshPoints);
  Pstream::gatherList(gatheredPoints);

  // Gather all faces
  List<faceList> gatheredFaces(Pstream::nProcs());
  gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
  forAll (gatheredFaces[Pstream::myProcNo()], i) {
    inplaceRenumber(pointToGlobal, gatheredFaces[Pstream::myProcNo()][i]);
  }
  Pstream::gatherList(gatheredFaces);

  // Gather all ZoneIDs
  List<labelList> gatheredZones(Pstream::nProcs());

  gatheredZones[Pstream::myProcNo()].transfer(compactZones);

  Pstream::gatherList(gatheredZones);

  // On master combine all points, faces, zones
  if (Pstream::master()) {
    pointField allPoints = ListListOps::combine<pointField>(
        gatheredPoints, accessOp<pointField>());
    gatheredPoints.clear();

    faceList allFaces =
        ListListOps::combine<faceList>(gatheredFaces, accessOp<faceList>());
    gatheredFaces.clear();

    labelList allZones =
        ListListOps::combine<labelList>(gatheredZones, accessOp<labelList>());
    gatheredZones.clear();

    // Zones
    surfZoneIdentifierList surfZones(compactZoneID.size());
    forAllConstIters(compactZoneID, iter) {
      surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
      Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
          << endl;
    }

    UnsortedMeshedSurface<face> unsortedFace(std::move(allPoints),
                                             std::move(allFaces),
                                             std::move(allZones),
                                             surfZones);

    MeshedSurface<face> sortedFace(unsortedFace);

    fileName globalCasePath(runTime.processorCase()
                                ? runTime.path() / ".." / outFileName
                                : runTime.path() / outFileName);
    globalCasePath.clean();

    Info << "Writing merged surface to " << globalCasePath << endl;

    sortedFace.write(globalCasePath);
  }
}

Foam::label snappymeshGen::checkAlignment(const Foam::polyMesh& mesh,
                                          const Foam::scalar tol,
                                          Foam::Ostream& os) {
  // Check all edges aligned with one of the coordinate axes
  const faceList &faces = mesh.faces();
  const pointField &points = mesh.points();

  label nUnaligned = 0;

  forAll(faces, facei) {
    const face &f = faces[facei];
    forAll(f, fp) {
      label fp1 = f.fcIndex(fp);
      const linePointRef e(edge(f[fp], f[fp1]).line(points));
      const vector v(e.vec());
      const scalar magV(mag(v));
      if (magV > ROOTVSMALL) {
        for (direction dir = 0; dir < pTraits<vector>::nComponents; ++dir) {
          const scalar s(mag(v[dir]));
          if (s > magV * tol && s < magV * (1 - tol)) {
            ++nUnaligned;
            break;
          }
        }
      }
    }
  }

  reduce(nUnaligned, sumOp<label>());

  if (nUnaligned) {
    os << "Initial mesh has " << nUnaligned
       << " edges unaligned with any of the coordinate axes" << nl << endl;
  }
  return nUnaligned;
}

int snappymeshGen::createMeshFromSTL(const char *fname) {
  using namespace Foam;

  int argc = 1;
  char **argv = new char *[2];
  argv[0] = new char[100];
  strcpy(argv[0], "NONE");
  Foam::argList args(argc, argv);
  Foam::Info << "Create time\n" << Foam::endl;
  Foam::argList::noParallel();

  Foam::fileName one = ".";
  Foam::fileName two = ".";
  Time runTime(Time::controlDictName, one, two);

  const bool overwrite = true;
  const bool dryRun = false;
  const bool checkGeometry = false;
  const bool surfaceSimplify = false;

  autoPtr<fvMesh> meshPtr;
  {
    Foam::Info << "Create mesh for time = " << runTime.timeName() << Foam::nl
               << Foam::endl;

    meshPtr.set(new fvMesh(Foam::IOobject(Foam::fvMesh::defaultRegion,
                                          runTime.timeName(), runTime,
                                          Foam::IOobject::MUST_READ)));
  }

  fvMesh &mesh = meshPtr();

  Info << "Read mesh in = " << runTime.cpuTimeIncrement() << " s" << endl;

  // Check patches and faceZones are synchronised
  mesh.boundaryMesh().checkParallelSync(true);
  meshRefinement::checkCoupledFaceZones(mesh);

  // Read meshing dictionary
  const word dictName("snappyHexMeshDict");
#include <setSystemMeshDictionaryIO.H>
  const IOdictionary meshDict(dictIO);

  // all surface geometry
  const dictionary& geometryDict =
      meshRefinement::subDict(meshDict, "geometry", dryRun);

  // refinement parameters
  const dictionary& refineDict =
      meshRefinement::subDict(meshDict, "castellatedMeshControls", dryRun);

  // mesh motion and mesh quality parameters
  const dictionary& motionDict =
      meshRefinement::subDict(meshDict, "meshQualityControls", dryRun);

  // snap-to-surface parameters
  const dictionary& snapDict =
      meshRefinement::subDict(meshDict, "snapControls", dryRun);

  // layer addition parameters
  const dictionary& layerDict =
      meshRefinement::subDict(meshDict, "addLayersControls", dryRun);

  // absolute merge distance
  const scalar mergeDist = getMergeDistance(
      mesh,
      meshRefinement::get<scalar>(meshDict,"mergeTolerance",dryRun),dryRun);

  const bool keepPatches(meshDict.getOrDefault("keepPatches", false));

  // format to be used for writing lines
  const word setFormat
  (
    meshDict.getOrDefault<word>
    (
      "setFormat",
      vtkSetWriter<scalar>::typeName
    )
  );
  const autoPtr<writer<scalar>> setFormatter
  (
    writer<scalar>::New(setFormat)
  );

  const scalar maxSizeRatio
  (
    meshDict.getOrDefault<scalar>("maxSizeRatio", 100)
  );

  // Read decomposePar dictionary
  dictionary decomposeDict;
  {
    if (Pstream::parRun()) {
      IOdictionary* dictPtr = new IOdictionary
      (
        IOobject::selectIO
        (
          IOobject
          (
            decompositionModel::canonicalName,
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
          ),
          args.getOrDefault<fileName>("decomposeParDict", "")
        )
      );

      // Store it on the object registry, but to be found it must also
      // have the expected "decomposeParDict" name.

      dictPtr->rename(decompositionModel::canonicalName);
      runTime.store(dictPtr);

      decomposeDict = *dictPtr;
    } else {
      decomposeDict.add("method", "none");
      decomposeDict.add("numberOfSubdomains", 1);
    }
  }

  // Debug
  // ~~~~~

  // Set debug level
  meshRefinement::debugType debugLevel =
      meshRefinement::debugType(meshDict.getOrDefault<label>("debug", 0));
  {
    wordList flags;
    if (meshDict.readIfPresent("debugFlags", flags)) {
      debugLevel = meshRefinement::debugType(
          meshRefinement::readFlags(meshRefinement::debugTypeNames, flags));
    }
  }
  if (debugLevel > 0) {
    meshRefinement::debug = debugLevel;
    snappyRefineDriver::debug = debugLevel;
    snappySnapDriver::debug = debugLevel;
    snappyLayerDriver::debug = debugLevel;
  }

  // Set file writing level
  {
    wordList flags;
    if (meshDict.readIfPresent("writeFlags", flags)) {
      meshRefinement::writeLevel(meshRefinement::writeType(
          meshRefinement::readFlags(meshRefinement::writeTypeNames, flags)));
    }
  }

  // for the impatient who want to see some output files:
  profiling::writeNow();

  // Read geometry
  // ~~~~~~~~~~~~~
  searchableSurfaces allGeometry(
      IOobject("abc",                   // dummy name
               mesh.time().constant(),  // instance
               // mesh.time().findInstance("triSurface", word::null),// instance
               "triSurface",  // local
               mesh.time(),   // registry
               IOobject::MUST_READ, IOobject::NO_WRITE),
      geometryDict, meshDict.getOrDefault("singleRegionName", true));

  // Read refinement surfaces
  // ~~~~~~~~~~~~~~~~~~~~~~~~

  autoPtr<refinementSurfaces> surfacesPtr;

  Info << "Reading refinement surfaces." << endl;
  
  surfacesPtr.reset(
    new refinementSurfaces(allGeometry,meshRefinement::subDict(
      refineDict,
      "refinementSurfaces",
      dryRun), refineDict.getOrDefault("gapLevelIncrement", 0), dryRun));

  Info << "Read refinement surfaces in = " << mesh.time().cpuTimeIncrement()
       << " s" << nl << endl;

  refinementSurfaces &surfaces = surfacesPtr();

  if (dryRun) {
    // Check geometry to mesh bounding box
    Info<< "Checking for geometry size relative to mesh." << endl;
    const boundBox& meshBb = mesh.bounds();
    forAll(allGeometry, geomi)
    {
      const searchableSurface& s = allGeometry[geomi];
      const boundBox& bb = s.bounds();

      scalar ratio = bb.mag() / meshBb.mag();
      if (ratio > maxSizeRatio || ratio < 1.0/maxSizeRatio)
      {
        Warning
            << "    " << allGeometry.names()[geomi]
            << " bounds differ from mesh"
            << " by more than a factor " << maxSizeRatio << ":" << nl
            << "        bounding box      : " << bb << nl
            << "        mesh bounding box : " << meshBb
            << endl;
      }
      if (!meshBb.contains(bb))
      {
        Warning
            << "    " << allGeometry.names()[geomi]
            << " bounds not fully contained in mesh" << nl
            << "        bounding box      : " << bb << nl
            << "        mesh bounding box : " << meshBb
            << endl;
      }
    }
    Info<< endl;
  }

  // Read refinement shells
  // ~~~~~~~~~~~~~~~~~~~~~~

  Info << "Reading refinement shells." << endl;
  shellSurfaces shells(
    allGeometry,
    meshRefinement::subDict(refineDict, "refinementRegions", dryRun),
    dryRun
  );
  Info << "Read refinement shells in = " << mesh.time().cpuTimeIncrement()
       << " s" << nl << endl;

  Info << "Setting refinement level of surface to be consistent"
       << " with shells." << endl;
  surfaces.setMinLevelFields(shells);
  Info << "Checked shell refinement in = " << mesh.time().cpuTimeIncrement()
       << " s" << nl << endl;

  // Optionally read limit shells
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  const dictionary limitDict(refineDict.subOrEmptyDict("limitRegions"));

  if (!limitDict.empty()) {
    Info << "Reading limit shells." << endl;
  }

  shellSurfaces limitShells(allGeometry, limitDict, dryRun);

  if (!limitDict.empty()) {
    Info << "Read limit shells in = " << mesh.time().cpuTimeIncrement() << " s"
         << nl << endl;
  }

  if (dryRun) {
    // Check for use of all geometry
    const wordList &allGeomNames = allGeometry.names();

    labelHashSet unusedGeometries(identity(allGeomNames.size()));
    unusedGeometries.erase(surfaces.surfaces());
    unusedGeometries.erase(shells.shells());
    unusedGeometries.erase(limitShells.shells());

    if (unusedGeometries.size()) {
      IOWarningInFunction(geometryDict)
          << "The following geometry entries are not used:" << nl;
      for (const label geomi : unusedGeometries) {
        Info << "    " << allGeomNames[geomi] << nl;
      }
      Info << endl;
    }
  }

  // Read feature meshes
  // ~~~~~~~~~~~~~~~~~~~

  Info << "Reading features." << endl;
  refinementFeatures features(mesh,
                              PtrList<dictionary>(meshRefinement::lookup(
                                  refineDict, "features", dryRun)),
                              dryRun);
  Info << "Read features in = " << mesh.time().cpuTimeIncrement() << " s" << nl
       << endl;

  if (dryRun) {
    // Check geometry to mesh bounding box
    Info << "Checking for line geometry size relative to surface geometry."
         << endl;

    OStringStream os;
    bool hasErrors =
        features.checkSizes(maxSizeRatio,  // const scalar maxRatio,
                            mesh.bounds(),
                            true,  // const bool report,
                            os     // FatalIOError
        );
    if (hasErrors) {
      Warning << os.str() << endl;
    }
  }

  // Refinement engine
  // ~~~~~~~~~~~~~~~~~

  Info << nl << "Determining initial surface intersections" << nl
       << "-----------------------------------------" << nl << endl;

  // Main refinement engine
  meshRefinement meshRefiner(
      mesh,
      mergeDist,    // tolerance used in sorting coordinates
      overwrite,    // overwrite mesh files?
      surfaces,     // for surface intersection refinement
      features,     // for feature edges/point based refinement
      shells,       // for volume (inside/outside) refinement
      limitShells,  // limit of volume refinement
      labelList(),  // initial faces to test
      dryRun);

  if (!dryRun) {
    meshRefiner.updateIntersections(identity(mesh.nFaces()));
    Info << "Calculated surface intersections in = "
         << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
  }

  // Some stats
  meshRefiner.printMeshInfo(debugLevel, "Initial mesh");

  meshRefiner.write(
      meshRefinement::debugType(debugLevel & meshRefinement::OBJINTERSECTIONS),
      meshRefinement::writeType(0),
      mesh.time().path() / meshRefiner.timeName());

  // Refinement parameters
  const refinementParameters refineParams(refineDict, dryRun);

  // Snap parameters
  const snapParameters snapParams(snapDict, dryRun);

  // Add all the cellZones and faceZones
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // 1. cellZones relating to surface (faceZones added later)

  const labelList namedSurfaces(
      surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones()));

  labelList surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh(
      surfaces.surfZones(), namedSurfaces, mesh);

  // 2. cellZones relating to locations

  refineParams.addCellZonesToMesh(mesh);

  // Add all the surface regions as patches
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //- Global surface region to patch (non faceZone surface) or patches
  //  (faceZone surfaces)
  labelList globalToMasterPatch;
  labelList globalToSlavePatch;

  {
    Info << nl << "Adding patches for surface regions" << nl
         << "----------------------------------" << nl << endl;

    // From global region number to mesh patch.
    globalToMasterPatch.setSize(surfaces.nRegions(), -1);
    globalToSlavePatch.setSize(surfaces.nRegions(), -1);

    if (!dryRun) {
      Info << setf(ios_base::left) << Foam::setw(6) << "Patch" << Foam::setw(20)
           << "Type" << Foam::setw(30) << "Region" << nl << Foam::setw(6)
           << "-----" << Foam::setw(20) << "----" << Foam::setw(30) << "------"
           << endl;
    }

    const labelList &surfaceGeometry = surfaces.surfaces();
    const PtrList<dictionary> &surfacePatchInfo = surfaces.patchInfo();
    const polyBoundaryMesh &pbm = mesh.boundaryMesh();

    forAll(surfaceGeometry, surfi) {
      label geomi = surfaceGeometry[surfi];

      const wordList &regNames = allGeometry.regionNames()[geomi];

      if (!dryRun) {
        Info << surfaces.names()[surfi] << ':' << nl << nl;
      }

      const wordList &fzNames = surfaces.surfZones()[surfi].faceZoneNames();

      if (fzNames.empty()) {
        // 'Normal' surface
        forAll(regNames, i) {
          label globalRegioni = surfaces.globalRegion(surfi, i);

          label patchi;

          if (surfacePatchInfo.set(globalRegioni)) {
            patchi = meshRefiner.addMeshedPatch(
                regNames[i], surfacePatchInfo[globalRegioni]);
          } else {
            dictionary patchInfo;
            patchInfo.set("type", wallPolyPatch::typeName);

            patchi = meshRefiner.addMeshedPatch(regNames[i], patchInfo);
          }

          if (!dryRun) {
            Info << setf(ios_base::left) << Foam::setw(6) << patchi
                 << Foam::setw(20) << pbm[patchi].type() << Foam::setw(30)
                 << regNames[i] << nl;
          }

          globalToMasterPatch[globalRegioni] = patchi;
          globalToSlavePatch[globalRegioni] = patchi;
        }
      } else {
        // Zoned surface
        forAll(regNames, i) {
          label globalRegioni = surfaces.globalRegion(surfi, i);

          // Add master side patch
          {
            label patchi;

            if (surfacePatchInfo.set(globalRegioni)) {
              patchi = meshRefiner.addMeshedPatch(
                  regNames[i], surfacePatchInfo[globalRegioni]);
            } else {
              dictionary patchInfo;
              patchInfo.set("type", wallPolyPatch::typeName);

              patchi = meshRefiner.addMeshedPatch(regNames[i], patchInfo);
            }

            if (!dryRun) {
              Info << setf(ios_base::left) << Foam::setw(6) << patchi
                   << Foam::setw(20) << pbm[patchi].type() << Foam::setw(30)
                   << regNames[i] << nl;
            }

            globalToMasterPatch[globalRegioni] = patchi;
          }
          // Add slave side patch
          {
            const word slaveName = regNames[i] + "_slave";
            label patchi;

            if (surfacePatchInfo.set(globalRegioni)) {
              patchi = meshRefiner.addMeshedPatch(
                  slaveName, surfacePatchInfo[globalRegioni]);
            } else {
              dictionary patchInfo;
              patchInfo.set("type", wallPolyPatch::typeName);

              patchi = meshRefiner.addMeshedPatch(slaveName, patchInfo);
            }

            if (!dryRun) {
              Info << setf(ios_base::left) << Foam::setw(6) << patchi
                   << Foam::setw(20) << pbm[patchi].type() << Foam::setw(30)
                   << slaveName << nl;
            }

            globalToSlavePatch[globalRegioni] = patchi;
          }
        }

        // For now: have single faceZone per surface. Use first
        // region in surface for patch for zoning
        if (regNames.size()) {
          forAll(fzNames, fzi) {
            const word &fzName = fzNames[fzi];
            label globalRegioni = surfaces.globalRegion(surfi, fzi);

            meshRefiner.addFaceZone(
                fzName, pbm[globalToMasterPatch[globalRegioni]].name(),
                pbm[globalToSlavePatch[globalRegioni]].name(),
                surfaces.surfZones()[surfi].faceType());
          }
        }
      }

      if (!dryRun) {
        Info << nl;
      }
    }
    Info << "Added patches in = " << mesh.time().cpuTimeIncrement() << " s"
         << nl << endl;
  }

  // Add all information for all the remaining faceZones
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  HashTable<Pair<word>> faceZoneToPatches;
  forAll(mesh.faceZones(), zonei) {
    const word &fzName = mesh.faceZones()[zonei].name();

    label mpI, spI;
    surfaceZonesInfo::faceZoneType fzType;
    bool hasInfo = meshRefiner.getFaceZoneInfo(fzName, mpI, spI, fzType);

    if (!hasInfo) {
      // faceZone does not originate from a surface but presumably
      // from a cellZone pair instead
      string::size_type i = fzName.find("_to_");
      if (i != string::npos) {
        word cz0 = fzName.substr(0, i);
        word cz1 = fzName.substr(i + 4, fzName.size() - i + 4);
        word slaveName(cz1 + "_to_" + cz0);
        faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
      } else {
        // Add as fzName + fzName_slave
        const word slaveName = fzName + "_slave";
        faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
      }
    }
  }

  if (faceZoneToPatches.size()) {
    snappyRefineDriver::addFaceZones(meshRefiner, refineParams,
                                     faceZoneToPatches);
  }

  // Re-do intersections on meshed boundaries since they use an extrapolated
  // other side
  {
    const labelList adaptPatchIDs(meshRefiner.meshedPatches());

    const polyBoundaryMesh &pbm = mesh.boundaryMesh();

    label nFaces = 0;
    forAll(adaptPatchIDs, i) { nFaces += pbm[adaptPatchIDs[i]].size(); }

    labelList faceLabels(nFaces);
    nFaces = 0;
    forAll(adaptPatchIDs, i) {
      const polyPatch &pp = pbm[adaptPatchIDs[i]];
      forAll(pp, i) { faceLabels[nFaces++] = pp.start() + i; }
    }
    meshRefiner.updateIntersections(faceLabels);
  }

  // Parallel
  // ~~~~~~~~

  // Construct decomposition engine. Note: cannot use decompositionModel
  // MeshObject since we're clearing out the mesh inside the mesh generation.
  autoPtr<decompositionMethod> decomposerPtr(
      decompositionMethod::New(decomposeDict));
  decompositionMethod &decomposer = *decomposerPtr;

  if (Pstream::parRun() && !decomposer.parallelAware()) {
    FatalErrorInFunction << "You have selected decomposition method "
                         << decomposer.typeName
                         << " which is not parallel aware." << endl
                         << "Please select one that is (hierarchical, ptscotch)"
                         << exit(FatalError);
  }

  // Mesh distribution engine (uses tolerance to reconstruct meshes)
  fvMeshDistribute distributor(mesh, mergeDist);

  // Now do the real work -refinement -snapping -layers
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  const bool wantRefine(
      meshRefinement::get<bool>(meshDict, "castellatedMesh", dryRun));
  const bool wantSnap(meshRefinement::get<bool>(meshDict, "snap", dryRun));
  const bool wantLayers(
      meshRefinement::get<bool>(meshDict, "addLayers", dryRun));

  if (dryRun) {
    string errorMsg(FatalError.message());
    string IOerrorMsg(FatalIOError.message());

    if (errorMsg.size() || IOerrorMsg.size()) {
      // errorMsg = "[dryRun] " + errorMsg;
      // errorMsg.replaceAll("\n", "\n[dryRun] ");
      // IOerrorMsg = "[dryRun] " + IOerrorMsg;
      // IOerrorMsg.replaceAll("\n", "\n[dryRun] ");

      Warning << nl << "Missing/incorrect required dictionary entries:" << nl
              << nl << IOerrorMsg.c_str() << nl << errorMsg.c_str() << nl << nl
              << "Exiting dry-run" << nl << endl;

      FatalError.clear();
      FatalIOError.clear();

      return 0;
    }
  }

  // How to treat co-planar faces
  meshRefinement::FaceMergeType mergeType =
      meshRefinement::FaceMergeType::GEOMETRIC;
  {
    const bool mergePatchFaces(meshDict.getOrDefault("mergePatchFaces", true));

    if (!mergePatchFaces) {
      Info << "Not merging patch-faces of cell to preserve"
           << " (split)hex cell shape." << nl << endl;
      mergeType = meshRefinement::FaceMergeType::NONE;
    } else {
      const bool mergeAcrossPatches(
          meshDict.getOrDefault("mergeAcrossPatches", false));

      if (mergeAcrossPatches) {
        Info << "Merging co-planar patch-faces of cells"
             << ", regardless of patch assignment" << nl << endl;
        mergeType = meshRefinement::FaceMergeType::IGNOREPATCH;
      }
    }
  }

  // Layer addition parameters
  const layerParameters layerParams(layerDict, mesh.boundaryMesh());

  if (wantRefine) {
    cpuTime timer;

    snappyRefineDriver refineDriver(meshRefiner, decomposer, distributor,
                                    globalToMasterPatch, globalToSlavePatch,
                                    setFormatter,dryRun);

    if (!overwrite && !debugLevel) {
      const_cast<Time &>(mesh.time())++;
    }

    refineDriver.doRefine(refineDict, refineParams, snapParams,
                          refineParams.handleSnapProblems(),mergeType,motionDict);

    if (!keepPatches && !wantSnap && !wantLayers) {
      fvMeshTools::removeEmptyPatches(mesh,true);
    }

    writeMesh("Refined mesh", meshRefiner, debugLevel,
              meshRefinement::writeLevel());

    Info << "Mesh refined in = " << timer.cpuTimeIncrement() << " s." << endl;

    profiling::writeNow();
  }

  if (wantSnap) {
    cpuTime timer;
    snappySnapDriver snapDriver(meshRefiner, globalToMasterPatch,
                                globalToSlavePatch,false);


    if (!overwrite && !debugLevel) {
      const_cast<Time &>(mesh.time())++;
    }

    // Use the resolveFeatureAngle from the refinement parameters
    scalar curvature = refineParams.curvature();
    scalar planarAngle = refineParams.planarAngle();

    snapDriver.doSnap(snapDict, motionDict, mergeType, curvature, planarAngle, snapParams);

    if (!keepPatches && !wantLayers) {
      fvMeshTools::removeEmptyPatches(mesh,true);
    }

    writeMesh("Snapped mesh", meshRefiner, debugLevel,
              meshRefinement::writeLevel());

    Info << "Mesh snapped in = " << timer.cpuTimeIncrement() << " s." << endl;

    profiling::writeNow();
  }

  if (wantLayers) {
    cpuTime timer;

    snappyLayerDriver layerDriver(meshRefiner, globalToMasterPatch,
                                  globalToSlavePatch,dryRun);

    // Use the maxLocalCells from the refinement parameters
    bool preBalance = returnReduce(
        (mesh.nCells() >= refineParams.maxLocalCells()), orOp<bool>());

    if (!overwrite && !debugLevel) {
      const_cast<Time &>(mesh.time())++;
    }

    layerDriver.doLayers(layerDict, motionDict, layerParams, mergeType, preBalance,
                         decomposer, distributor);

    if (!keepPatches) {
      fvMeshTools::removeEmptyPatches(mesh,true);
    }

    writeMesh("Layer mesh", meshRefiner, debugLevel,
              meshRefinement::writeLevel());

    Info << "Layers added in = " << timer.cpuTimeIncrement() << " s." << endl;
    profiling::writeNow();
  }

  {
    addProfiling(checkMesh, "snappyHexMesh::checkMesh");
    // Check final mesh
    Info << "Checking final mesh ..." << endl;
    faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces() / 100);

    motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces, dryRun);
    const label nErrors = returnReduce(wrongFaces.size(), sumOp<label>());

    if (nErrors > 0) {
      Info << "Finished meshing with " << nErrors << " illegal faces"
           << " (concave, zero area or negative cell pyramid volume)" << endl;
      wrongFaces.write();
    } else {
      Info << "Finished meshing without any errors" << endl;
    }
    profiling::writeNow();
  }

    profiling::writeNow();

  Info << "Finished meshing in = " << runTime.elapsedCpuTime() << " s." << endl;

  Info << "End\n" << endl;

  readSnappyFoamMesh();
}

void snappymeshGen::readSnappyFoamMesh() {
  Foam::Info << "Create time\n" << Foam::endl;
  Foam::argList::noParallel();

  Foam::fileName one = ".";
  Foam::fileName two = ".";
  Time runTime(Time::controlDictName, one, two);
  // reading mesh database and converting
  Foam::word regionName;
  regionName = Foam::fvMesh::defaultRegion;
  Foam::Info << "Create mesh for time = " << runTime.timeName() << Foam::nl
             << Foam::endl;

  _fmesh = new Foam::fvMesh(Foam::IOobject(regionName, runTime.timeName(),
                                           runTime, Foam::IOobject::MUST_READ));

  genMshDB();
}

void snappymeshGen::genMshDB() {
  // Obtaining the mesh data and generating requested output file
  std::cout << "Number of points " << _fmesh->nPoints() << std::endl;
  std::cout << "Number of cells " << _fmesh->nCells() << std::endl;

  // declare vtk dataset
  vtkSmartPointer<vtkUnstructuredGrid> dataSet_tmp =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

  // decomposition
  Foam::foamVTKTopo::decomposePoly = false;

  // creating equivalent vtk topology from fvMesh
  // by default polyhedral cells will be decomposed to
  // tets and pyramids. Additional points will be added
  // to underlying fvMesh.
  std::cout << "Performing topological decomposition.\n";
  Foam::foamVTKTopo topo(*_fmesh);

  // point coordinates
  Foam::pointField pf = _fmesh->points();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (int ipt = 0; ipt < _fmesh->nPoints(); ipt++)
    points->InsertNextPoint(pf[ipt].x(), pf[ipt].y(), pf[ipt].z());
  dataSet_tmp->SetPoints(points);

  // cell types
  std::vector<int> pntIds;
  int nCelPnts = 0;
  for (int icl = 0; icl < topo.vertLabels().size(); icl++) {
    if (topo.cellTypes()[icl] != VTK_POLYHEDRON) {
      nCelPnts = topo.vertLabels()[icl].size();
      pntIds.resize(nCelPnts, -1);
      for (int ip = 0; ip < nCelPnts; ip++)
        pntIds[ip] = topo.vertLabels()[icl][ip];
      createVtkCell(dataSet_tmp, topo.cellTypes()[icl], pntIds);
    } else {
      // polyhedral cells treated differently in vtk
      // faces should be defined for them
      int nFace = topo.vertLabels()[icl][0];
      vtkSmartPointer<vtkCellArray> faces =
          vtkSmartPointer<vtkCellArray>::New();
      std::vector<vtkIdType> faceIds;
      std::set<vtkIdType> pntIds;
      int dataId = 1;
      for (int iFace = 0; iFace < nFace; iFace++) {
        faceIds.clear();
        int nFaceId = topo.vertLabels()[icl][dataId++];
        int pntId;
        for (int ifid = 0; ifid < nFaceId; ifid++) {
          pntId = topo.vertLabels()[icl][dataId++];
          pntIds.insert(pntId);
          faceIds.push_back(pntId);
        }
        faces->InsertNextCell(nFaceId, &faceIds[0]);
      }
      std::vector<vtkIdType> pntIdsVec(pntIds.begin(), pntIds.end());
      dataSet_tmp->InsertNextCell(VTK_POLYHEDRON, pntIdsVec.size(),
                                  &pntIdsVec[0], nFace, faces->GetPointer());
    }
  }
  dataSet = dataSet_tmp;
}

void snappymeshGen::createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                                  const int cellType,
                                  std::vector<int> &vrtIds) {
  vtkSmartPointer<vtkIdList> vtkCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkCellIds->SetNumberOfIds(vrtIds.size());
  for (auto pit = vrtIds.begin(); pit != vrtIds.end(); pit++)
    vtkCellIds->SetId(pit - vrtIds.begin(), *pit);
  dataSet->InsertNextCell(cellType, vtkCellIds);
}
