#include <boost/filesystem.hpp>
#include <iostream>
#include <set>
#include <string>

#include "Mesh/foamGeoMesh.H"
#include "Mesh/geoMeshBase.H"
#include "Mesh/geoMeshFactory.H"
#include "MeshGeneration/snappymeshGen.H"
#include "MeshGeneration/snappymeshParams.H"

// openfoam headers
#include <getDicts.H>
#include <snappyMesh.H>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// This class has all the methods from snappyHexMesh. Currently, Surface
// Simplify method is not enabled.

snappymeshGen::snappymeshGen()  // default constructor definition
{
  // Default Meshing Parameters
  params_ = new snappymeshParams();
  defaults = true;

  // Initialization Tasks
  initialize();
}

// constructor definition with parameters input
snappymeshGen::snappymeshGen(snappymeshParams *params)
    : defaults(false), params_(params) {
  initialize();  // Foam Initialization
}

snappymeshGen::~snappymeshGen() {}

void snappymeshGen::initialize() {
  // create dictionaries needed
  bool writeDicts;
  if (params_->isPackMesh)
    writeDicts = true;
  else
    writeDicts = false;

  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  controlDict_ = initFoam->createControlDict(writeDicts);
  fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  fvSolution_ = initFoam->createFvSolution(writeDicts);

  // create dictionaries needed in memory
  createSnappyDict(writeDicts);

  runTime_ =
      std::unique_ptr<Foam::Time>(new Foam::Time(controlDict_.get(), ".", "."));

  //- 2d cartesian mesher cannot be run in parallel
  Foam::argList::noParallel();
}

int snappymeshGen::createMeshFromSTL(const char *fname) {
  auto objMsh = snappyMesh();

  bool writeMsh;
  if (params_->isPackMesh)
    writeMsh = true;
  else
    writeMsh = false;

  fmesh_ = std::unique_ptr<Foam::fvMesh>(
      objMsh.generate(runTime_, snappyMshDict_, nullptr, writeMsh));

  // Create foamGeoMesh
  if (writeMsh) {
    gmData.reset();
    gmData = std::unique_ptr<NEM::MSH::geoMeshBase>(NEM::MSH::Read(".foam"));
  } else {
    auto fgm_ = std::unique_ptr<NEM::MSH::foamGeoMesh>(
        new NEM::MSH::foamGeoMesh(fmesh_.get(), ""));
    gmData = std::unique_ptr<NEM::MSH::geoMeshBase>(
        dynamic_cast<NEM::MSH::geoMeshBase *>(fgm_.get()));
  }

  return 0;
}

void snappymeshGen::createSnappyDict(const bool &write) {
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
  if ((params_->withCastMesh)) {
    contText += "\ncastellatedMesh  true;\n";
  } else {
    contText += "\ncastellatedMesh  false;\n";
  }

  if ((params_->withSnap)) {
    contText += "\nsnap       true;\n";
  } else {
    contText += "\nsnap       false;\n";
  }

  if ((params_->withLayers)) {
    contText += "\naddLayers      true;\n";
  } else {
    contText += "\naddLayers      false;\n";
  }

  // Geometry declaration
  contText += "\ngeometry\n{\n";
  contText += "\t" + (params_->geomFileName);
  contText += "\n\t{\n\t\ttype\ttriSurfaceMesh;\n";

  if (params_->geoDef.withMultiPatches) {
    if (params_->geoDef.stlPatchDefs.empty()) {
      std::cerr << "Error reading in custom patches from JSON input!"
                << std::endl;
      throw;
    }
    contText += "\t\tregions\n\t\t{\n";
    for (auto pt = (params_->geoDef.stlPatchDefs).begin();
         pt != (params_->geoDef.stlPatchDefs).end(); pt++) {
      contText += "\t\t\t" + pt->STLPatchName + "  \t{name " +
                  pt->snappyPatchName + ";}\n";
    }
    contText += "\t\t}\n";
  } else {
    contText =
        contText + "\t\tname\t" + (params_->geoDef.singleSolidPatch) + ";\n";
  }

  // Searchable shapes go here
  if (!params_->geoDef.srchShape.empty()) {
    for (const auto &shape : params_->geoDef.srchShape) {
      contText += "\t\t" + (shape->patchName);
      contText += "\n\t\t{";
      contText += "\n\t\t\ttype " + (shape->getType()) + ";\n";

      if (auto box = std::dynamic_pointer_cast<shmSearchableBox>(shape)) {
        contText += "\t\t\tmin (" + std::to_string(box->minBound.at(0)) + " " +
                    std::to_string(box->minBound.at(1)) + " " +
                    std::to_string(box->minBound.at(2)) + ");\n";
        contText += "\t\t\tmax (" + std::to_string(box->maxBound.at(0)) + " " +
                    std::to_string(box->maxBound.at(1)) + " " +
                    std::to_string(box->maxBound.at(2)) + ");\n";
      } else if (auto sphere =
                     std::dynamic_pointer_cast<shmSearchableSphere>(shape)) {
        contText += "\t\t\tcentre (" + std::to_string(sphere->center.at(0)) +
                    " " + std::to_string(sphere->center.at(1)) + " " +
                    std::to_string(sphere->center.at(2)) + ");\n";
        contText += "\t\t\tradius " + std::to_string(sphere->radius) + ";\n";
      } else if (auto cyl =
                     std::dynamic_pointer_cast<shmSearchableCylinder>(shape)) {
        contText += "\t\t\tpoint1 (" + std::to_string(cyl->axisPoint1.at(0)) +
                    " " + std::to_string(cyl->axisPoint1.at(1)) + " " +
                    std::to_string(cyl->axisPoint1.at(2)) + ");\n";
        contText += "\t\t\tpoint2 (" + std::to_string(cyl->axisPoint2.at(0)) +
                    " " + std::to_string(cyl->axisPoint2.at(1)) + " " +
                    std::to_string(cyl->axisPoint2.at(2)) + ");\n";
        contText += "\t\t\tradius " + std::to_string(cyl->radius) + ";\n";
      }
      contText += "\n\t\t}\n";
    }
  }
  contText += "\t}\n";
  contText += "};\n";
  // End of geometry definition

  // Castellated Mesh Controls
  contText += "\n\ncastellatedMeshControls\n{\n";
  contText += "\tmaxLocalCells\t" +
              std::to_string(params_->castMeshControls.maxLCells) + ";\n";
  contText += "\tmaxGlobalCells\t" +
              std::to_string(params_->castMeshControls.maxGCells) + ";\n";
  contText += "\tnCellsBetweenLevels\t" +
              std::to_string(params_->castMeshControls.cellsBetnLvls) + ";\n";
  contText += "\tminRefinementCells\t" +
              std::to_string(params_->castMeshControls.minRefCells) + ";\n";

  // Features using emesh file
  contText += "\n\tfeatures\n\t(\n";
  if (!params_->castMeshControls.ftrEdge.empty()) {
    for (auto pt = (params_->castMeshControls.ftrEdge).begin();
         pt != (params_->castMeshControls.ftrEdge).end(); pt++) {
      contText += "\t\t{\n";
      contText += "\t\t\tfile \"" + pt->fileName + "\";\n";
      contText += "\t\t\tlevels ((" + std::to_string(pt->minLvl) + " " +
                  std::to_string(pt->maxLvl) + "));\n";
      contText += "\t\t}\n";
    }
  }
  contText += "\t);\n";

  contText += "\n\trefinementSurfaces\n\t{\n";
  if (params_->geoDef.withMultiPatches)
    contText += "\t\t" + (params_->geomFileName) + "\n\t\t{";
  else
    contText =
        contText + "\t\t" + (params_->geoDef.singleSolidPatch) + "\n\t\t{";
  contText += "\n\t\t\tlevel (" +
              std::to_string(params_->castMeshControls.refSurfLvlMin) + " " +
              std::to_string(params_->castMeshControls.refSurfLvlMax) + ");\n";

  if ((params_->castMeshControls.withCellZones)) {
    contText += "\t\t\tfaceZone\t" + (params_->geoDef.singleSolidPatch) + ";\n";
    contText += "\t\t\tcellZone\t" + (params_->geoDef.singleSolidPatch) + ";\n";
    contText += "\t\t\tcellZoneInside\tinside;\n";
  }

  if ((params_->geoDef.withMultiPatches)) {
    if ((!params_->castMeshControls.surfRefs.empty())) {
      contText += "\n\t\t\tregions\n";
      contText += "\t\t\t{";

      for (auto pt = (params_->castMeshControls.surfRefs).begin();
           pt != (params_->castMeshControls.surfRefs).end(); pt++) {
        if (pt->patchType == "NO") {
          contText += "\n\t\t\t\t" + (pt->refPatchName) + "\t{level (" +
                      std::to_string((int)pt->minLvl) + " " +
                      std::to_string(pt->maxLvl) +
                      "); patchInfo { type patch; }}";
        } else {
          contText += "\n\t\t\t\t" + (pt->refPatchName) + "\t{level (" +
                      std::to_string((int)pt->minLvl) + " " +
                      std::to_string(pt->maxLvl) + "); patchInfo { type " +
                      pt->patchType + "; }}";
        }
      }

      contText += "\n\t\t\t}\n";
    }
  }

  contText += "\t\t}\n";

  contText += "\t\tgapLevelIncrement " +
              std::to_string(params_->castMeshControls.castMeshGpLvl) + +";\n";

  contText += "\t}\n";

  contText += "\tresolveFeatureAngle\t" +
              std::to_string(params_->castMeshControls.featAngle) + ";\n";

  contText += "\tgapLevelIncrement\t" +
              std::to_string(params_->castMeshControls.gapLvlInc) + ";\n";

  contText += "\n\tplanarAngle " +
              std::to_string(params_->castMeshControls.planarAngle) + ";\n";

  contText += "\n\trefinementRegions\n\t{\n";

  if ((!params_->castMeshControls.geomRefs.empty())) {
    for (auto pt = (params_->castMeshControls.geomRefs).begin();
         pt != (params_->castMeshControls.geomRefs).end(); pt++) {
      contText += "\n\t\t" + (pt->patchName);
      contText += "\n\t\t{";
      contText += "\n\t\t\tmode " + (pt->mode) + ";\n";
      contText += "\t\t\tlevels ((" + std::to_string(pt->minLvl) + " " +
                  std::to_string(pt->maxLvl) + "));\n";
      contText += "\t\t}\n";
    }
  }

  contText += "\t}";

  contText += "\n\tlocationInMesh\t(" +
              std::to_string(params_->castMeshControls.locMesh.at(0)) + " " +
              std::to_string(params_->castMeshControls.locMesh.at(1)) + " " +
              std::to_string(params_->castMeshControls.locMesh.at(2)) + ");\n";
  if ((params_->castMeshControls.alwFreeZone) == 1)
    contText += "\tallowFreeStandingZoneFaces\ttrue;\n";
  if ((params_->castMeshControls.alwFreeZone) == 0)
    contText += "\tallowFreeStandingZoneFaces\tfalse;\n";
  contText += "}\n";

  // Snap Controls
  contText += "\n\nsnapControls\n{\n";
  contText += "\tnSmoothPatch\t" +
              std::to_string(params_->snapControls.nSmoothPatch) + ";\n";
  contText +=
      "\ttolerance\t" + std::to_string(params_->snapControls.tolerance) + ";\n";
  contText += "\tnSolveIter\t" +
              std::to_string(params_->snapControls.nSolveIter) + ";\n";
  contText += "\tnRelaxIter\t" +
              std::to_string(params_->snapControls.nRelaxIter) + ";\n";

  // New Features
  contText += "\tnFeatureSnapIter\t" +
              std::to_string(params_->snapControls.nFeatureSnapIter) + ";\n";

  if (params_->snapControls.implicitFeatureSnap)
    contText += "\timplicitFeatureSnap\t true;\n";
  else
    contText += "\timplicitFeatureSnap\t false;\n";

  if (params_->snapControls.explicitFeatureSnap)
    contText += "\texplicitFeatureSnap\t true;\n";
  else
    contText += "\texplicitFeatureSnap\t false;\n";

  if (params_->snapControls.multiRegionFeatureSnap)
    contText += "\tmultiRegionFeatureSnap\t true;\n";
  else
    contText += "\tmultiRegionFeatureSnap\t false;\n";

  contText += "}\n";

  // Layer Controls
  contText += "\n\naddLayersControls\n{\n";
  if ((params_->layerControls.relSize) == 1)
    contText += "\trelativeSizes\ttrue;\n";
  if ((params_->layerControls.relSize) == 0)
    contText += "\trelativeSizes\tfalse;\n";

  // Layers
  contText += "\tlayers\n\t{\n";
  if (!params_->layerControls.layerVec.empty()) {
    for (auto pt = (params_->layerControls.layerVec).begin();
         pt != (params_->layerControls.layerVec).end(); pt++) {
      contText += "\t\t" + pt->patchName + "\n\t\t{\n";
      contText +=
          "\t\t\tnSurfaceLayers " + std::to_string(pt->nSurfaceLayers) + ";\n";
      contText +=
          "\t\t\texpansionRatio " + std::to_string(pt->expansionRatio) + ";\n";
      contText += "\t\t\tfinalLayerThickness " +
                  std::to_string(pt->finalLayerThickness) + ";\n";
      if (pt->firstLyrThickness > 0)
        contText += "\t\t\tfirstLayerThickness " +
                    std::to_string(pt->firstLyrThickness) + ";\n";
      if (pt->thickness > 0)
        contText += "\t\t\tthickness " + std::to_string(pt->thickness) + ";\n";
      contText +=
          "\t\t\tminThickness " + std::to_string(pt->minThickness) + ";\n";
      contText += "\t\t}\n";
    }
  }
  contText += "\t}\n";

  contText += "\texpansionRatio\t" +
              std::to_string(params_->layerControls.expRatio) + ";\n";
  contText += "\tfinalLayerThickness\t" +
              std::to_string(params_->layerControls.finLThick) + ";\n";
  contText += "\tminThickness\t" +
              std::to_string(params_->layerControls.minThick) + ";\n";
  if (params_->layerControls.firstLyrThickness > 0)
    contText += "\tfirstLayerThickness\t" +
                std::to_string(params_->layerControls.firstLyrThickness) +
                ";\n";
  if (params_->layerControls.thickness > 0)
    contText += "\tthickness\t" +
                std::to_string(params_->layerControls.thickness) + ";\n";
  contText +=
      "\tnGrow\t" + std::to_string(params_->layerControls.nGrow) + ";\n";
  contText += "\tfeatureAngle\t" +
              std::to_string(params_->layerControls.featAngle) + ";\n";
  contText += "\tslipFeatureAngle\t" +
              std::to_string(params_->layerControls.slipFeatureAngle) + ";\n";
  contText += "\tnRelaxIter\t" +
              std::to_string(params_->layerControls.relaxIter) + ";\n";
  contText += "\tnSmoothSurfaceNormals\t" +
              std::to_string(params_->layerControls.smthSurfNorm) + ";\n";
  contText += "\tnSmoothNormals\t" +
              std::to_string(params_->layerControls.smthNorm) + ";\n";
  contText += "\tnSmoothThickness\t" +
              std::to_string(params_->layerControls.smthThick) + ";\n";
  contText += "\tmaxFaceThicknessRatio\t" +
              std::to_string(params_->layerControls.maxFcTR) + ";\n";
  contText += "\tmaxThicknessToMedialRatio\t" +
              std::to_string(params_->layerControls.maxThickTMR) + ";\n";
  contText += "\tminMedialAxisAngle\t" +
              std::to_string(params_->layerControls.minMedAngl) + ";\n";
  contText += "\tnBufferCellsNoExtrude\t" +
              std::to_string(params_->layerControls.bufferCells) + ";\n";
  contText +=
      "\tnLayerIter\t" + std::to_string(params_->layerControls.nIter) + ";\n";
  contText += "\tnRelaxedIter\t" +
              std::to_string(params_->layerControls.nRelaxedIter) + ";\n";
  if (params_->layerControls.nMedialAxisIter > 0)
    contText += "\tnMedialAxisIter\t" +
                std::to_string(params_->layerControls.nMedialAxisIter) + ";\n";
  if (params_->layerControls.nSmoothDisplacement > 0)
    contText += "\tnSmoothDisplacement\t" +
                std::to_string(params_->layerControls.nSmoothDisplacement) +
                ";\n";
  contText += "}\n";

  // Mesh Quality Controls
  std::ostringstream minimumVol;
  minimumVol << params_->qualityControls.minVol;
  std::ostringstream minimumTetQuality;
  minimumTetQuality << params_->qualityControls.minTetQ;

  contText += "\n\nmeshQualityControls\n{\n";
  contText += "\tmaxNonOrtho\t" +
              std::to_string(params_->qualityControls.maxNonOrtho) + ";\n";
  contText += "\tmaxBoundarySkewness\t" +
              std::to_string(params_->qualityControls.maxBndrySkew) + ";\n";
  contText += "\tmaxInternalSkewness\t" +
              std::to_string(params_->qualityControls.maxIntSkew) + ";\n";
  contText += "\tmaxConcave\t" +
              std::to_string(params_->qualityControls.maxConc) + ";\n";
  contText += "\tminVol\t" + minimumVol.str() + ";\n";
  contText += "\tminTetQuality\t" + minimumTetQuality.str() + ";\n";
  contText +=
      "\tminArea\t" + std::to_string(params_->qualityControls.minArea) + ";\n";
  contText += "\tminTwist\t" +
              std::to_string(params_->qualityControls.minTwist) + ";\n";
  contText += "\tminFaceWeight\t" +
              std::to_string(params_->qualityControls.minFaceW) + ";\n";
  contText += "\tminVolRatio\t" +
              std::to_string(params_->qualityControls.minVolRto) + ";\n";
  contText += "\tminDeterminant\t" +
              std::to_string(params_->qualityControls.minDet) + ";\n";
  contText += "\tminTriangleTwist\t" +
              std::to_string(params_->qualityControls.minTriTwist) + ";\n";
  contText += "\tnSmoothScale\t" +
              std::to_string(params_->qualityControls.smoothScale) + ";\n";
  contText += "\terrorReduction\t" +
              std::to_string(params_->qualityControls.errReduction) + ";\n";
  contText += "\trelaxed\n\t{\n";
  contText += "\t\tmaxNonOrtho \t " +
              std::to_string(params_->qualityControls.maxNonOrtho) + ";\n\t}\n";
  contText += "}\n";

  contText +=
      "\nmergeTolerance\t" + std::to_string(params_->mergeTol) + ";\n\n";

  contText += "writeFlags\n(\n";
  contText += "\tscalarLevels\n";
  contText += "\tlayerSets\n";
  contText += "\tlayerFields\n);\n";

  contText =
      contText +
      "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";

  // Keep this dictionary in memory
  Foam::dictionary tmptmpDuc =
      new Foam::dictionary(Foam::IStringStream(contText)(), true);
  snappyMshDict_ = std::unique_ptr<Foam::dictionary>(
      new Foam::dictionary("snappyHexMeshDict"));
  snappyMshDict_->merge(tmptmpDuc);

  // Write in traditional way due to formatting issues
  if (write) {
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
    contDict << contText;
    contDict.close();
  }
}
