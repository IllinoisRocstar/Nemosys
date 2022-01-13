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
#include <iostream>
#include <string>
#include "MeshGeneration/cfmeshGen.H"
#include "MeshGeneration/cfmeshParams.H"
#include "Mesh/geoMeshFactory.H"
#include "Mesh/geoMeshBase.H"
#include "Mesh/foamGeoMesh.H"

// vtk
#include <vtkIdList.h>

// openfoam headers
#include <fileName.H>
#include <getDicts.H>

// cfmesh headers
#include <cartesian2DMeshGenerator.H>
#include <triSurfaceDetectFeatureEdges.H>
#include <triSurfacePatchManipulator.H>
#include <triSurf.H>
#include <tetMeshGenerator.H>
#include <polyMeshGenModifier.H>
#include <meshOptimizer.H>
#include <cartesianMeshGenerator.H>
#include <voronoiMeshGenerator.H>

cfmeshGen::cfmeshGen() {
  // default meshing parameters
  params_ = new cfmeshParams();
  defaults = true;

  // Initialization tasks
  initialize();
}

cfmeshGen::cfmeshGen(cfmeshParams *params) : defaults(false), params_(params) {
  // Initialization tasks
  initialize();
}

cfmeshGen::~cfmeshGen() {}

void cfmeshGen::initialize() {
  // surface feature edge treatment
  if (params_->srfEdge.has_value())
    if (surfaceFeatureEdgeDetect()) {
      std::cerr << "A problem occured during edge detection step!\n";
      throw;
    }

  // create dictionaries needed in memory
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
  createMshDict(writeDicts);

  runTime_ =
      std::unique_ptr<Foam::Time>(new Foam::Time(controlDict_.get(), ".", "."));

  //- 2d cartesian mesher cannot be run in parallel
  Foam::argList::noParallel();
}

int cfmeshGen::createMeshFromSTL(const char *fname) {
  bool writeMsh;
  if (params_->isPackMesh)
    writeMsh = true;
  else
    writeMsh = false;

  // mesh generation and I/O
  Foam::Info << "Generating mesh with cfMesh engine" << Foam::endl;
  if (params_->generator == "cartesian2D") {
    Foam::Module::cartesian2DMeshGenerator cmg(*runTime_);
    std::cout << "ExecutionTime = " << runTime_->elapsedCpuTime() << " s\n"
              << "ClockTime = " << runTime_->elapsedClockTime() << " s"
              << std::endl;
    cmg.writeMesh();
  } else if (params_->generator == "tetMesh") {
    Foam::Module::tetMeshGenerator tmg(*runTime_);
    std::cout << "ExecutionTime = " << runTime_->elapsedCpuTime() << " s\n"
              << "ClockTime = " << runTime_->elapsedClockTime() << " s"
              << std::endl;
    tmg.writeMesh();

    // post-processing steps
    if (params_->improveMeshQuality.has_value()) improveMeshQuality();
  } else if (params_->generator == "cartesian3D") {
    Foam::Module::cartesianMeshGenerator cmg(*runTime_);
    std::cout << "ExecutionTime = " << runTime_->elapsedCpuTime() << " s\n"
              << "ClockTime = " << runTime_->elapsedClockTime() << " s"
              << std::endl;
    cmg.writeMesh();
  } else if (params_->generator == "polyMesh") {
    Foam::Module::voronoiMeshGenerator pmg(*runTime_);
    std::cout << "ExecutionTime = " << runTime_->elapsedCpuTime() << " s\n"
              << "ClockTime = " << runTime_->elapsedClockTime() << " s"
              << std::endl;
    pmg.writeMesh();
  } else {
    std::cerr << (params_->generator)
              << " is not a supported mesh generator.\n";
    throw;
  }

  // Create foamGeoMesh
  if (writeMsh) {
    gmData.reset();
    gmData = std::unique_ptr<NEM::MSH::geoMeshBase>(NEM::MSH::Read(".foam"));
  } else {
    auto fgm_ = std::unique_ptr<NEM::MSH::foamGeoMesh>(
        new NEM::MSH::foamGeoMesh(fmesh_.get(), ""));
    gmData = std::unique_ptr<NEM::MSH::geoMeshBase>(
        dynamic_cast<NEM::MSH::geoMeshBase *>(fgm_.get()));
    return 0;
  }

  return 0;
}

int cfmeshGen::surfaceFeatureEdgeDetect() {
  std::cout << "Performing surface feature edge detection.\n";
  std::string of = "./" + caseName + "_feature.ftr";
  Foam::fileName inFileName(params_->geomFilePath);
  Foam::fileName outFileName(of);

  if (outFileName == inFileName) {
    std::cerr << "Output file " << outFileName
              << " would overwrite the input file.\n";
    throw;
  }

  double tol = params_->srfEdge.value().srfEdgAng;
  std::cout << "Using " << tol << " deg angle\n";

  Foam::Module::triSurf originalSurface(inFileName);

  Foam::Module::triSurfaceDetectFeatureEdges edgeDetector(originalSurface, tol);
  edgeDetector.detectFeatureEdges();

  if (outFileName.ext() == "fms" || outFileName.ext() == "FMS") {
    std::cout << "Writing : " << outFileName << std::endl;
    originalSurface.writeSurface(outFileName);
  } else {
    Foam::Module::triSurfacePatchManipulator manipulator(originalSurface);
    const Foam::Module::triSurf *newSurfPtr = manipulator.surfaceWithPatches();

    std::cout << "Writing : " << outFileName << std::endl;
    newSurfPtr->writeSurface(outFileName);
    delete newSurfPtr;
  }

  // change cad input file
  params_->geomFilePath = of;

  return 0;
}

int cfmeshGen::improveMeshQuality() {
  std::cout << "Performing mesh quality improvements.\n";

  //- load the mesh from disk
  Foam::Module::polyMeshGen pmg(*runTime_);
  pmg.read();

  //- construct the smoother
  Foam::Module::meshOptimizer mOpt(pmg);

  const auto &meshQual = params_->improveMeshQuality.value();

  if ((meshQual.qltConCelSet) != "none") {
    //- lock cells in constrainedCellSet
    mOpt.lockCellsInSubset((meshQual.qltConCelSet));

    //- find boundary faces which shall be locked
    Foam::Module::labelLongList lockedBndFaces, selectedCells;

    const Foam::label sId = pmg.cellSubsetIndex((meshQual.qltConCelSet));
    pmg.cellsInSubset(sId, selectedCells);

    Foam::boolList activeCell(pmg.cells().size(), false);
    for (int iCl = 0; iCl < selectedCells.size(); iCl++)
      activeCell[selectedCells[iCl]] = true;
  }

  //- clear geometry information before volume smoothing
  pmg.clearAddressingData();

  //- perform optimisation using the laplace smoother and
  mOpt.optimizeMeshFV((meshQual.qltNLop), (meshQual.qltNLop),
                      (meshQual.qltNItr), (meshQual.qltNSrfItr));

  //- perform optimisation of worst quality faces
  mOpt.optimizeMeshFVBestQuality((meshQual.qltNLop), (meshQual.qltQltThr));

  //- check the mesh again and untangl bad regions if any of them exist
  mOpt.untangleMeshFV((meshQual.qltNLop), (meshQual.qltNItr),
                      (meshQual.qltNSrfItr));

  std::cout << "Finished optimization cycle\n";
  pmg.write();

  return 0;
}

void cfmeshGen::createMshDict(const bool &write) {
  meshDict_ =
      std::unique_ptr<Foam::dictionary>(new Foam::dictionary("meshDict"));

  Foam::dictionary fmfle("FoamFile");
  fmfle.add("version", "2.0");
  fmfle.add("format", "ascii");
  fmfle.add("class", "dictionary");
  fmfle.add("location", "\"system\"");
  fmfle.add("object", "meshDict");
  meshDict_->add("FoamFile", fmfle);

  if ((params_->maxCellSize) > 0)
    meshDict_->add("maxCellSize", params_->maxCellSize);

  if ((params_->minCellSize) > 0)
    meshDict_->add("minCellSize", params_->minCellSize);

  if ((params_->bndryCellSize) > 0)
    meshDict_->add("boundaryCellSize", params_->bndryCellSize);

  if ((params_->bndryCellSizeRefThk) > 0)
    meshDict_->add("boundaryCellSizeRefinementThickness",
                   params_->bndryCellSizeRefThk);

  if (params_->keepCellIB) meshDict_->add("keepCellsIntersectingBoundary", "1");

  if (params_->chkGluMsh) meshDict_->add("checkForGluedMesh", "1");

  if ((params_->alwDiscDomains)) meshDict_->add("allowDisconnectedDomains", 1);

  meshDict_->add("surfaceFile", params_->geomFilePath);

  // boundary layer
  if (params_->boundaryLayers.has_value()) {
    const auto &boundaryLayer = params_->boundaryLayers.value();
    Foam::dictionary boundaryLayrs("boundaryLayers");
    boundaryLayrs.add("nLayers", boundaryLayer.blNLyr);
    boundaryLayrs.add("thicknessRatio", boundaryLayer.blThkRto);
    if (boundaryLayer.maxFrstLyrThk > 0.0)
      boundaryLayrs.add("maxFirstLayerThickness", boundaryLayer.maxFrstLyrThk);

    if (boundaryLayer.alwDiscont) boundaryLayrs.add("allowDiscontinuity", 1);

    // boundary layer patches
    if (!boundaryLayer.blPatches.empty()) {
      Foam::dictionary patchBoundaryLayers("patchBoundaryLayers");
      for (auto pt = (boundaryLayer.blPatches).begin();
           pt != (boundaryLayer.blPatches).end(); pt++) {
        Foam::dictionary tmpDict1(pt->patchName);
        if (pt->alwDiscont) tmpDict1.add("allowDiscontinuity", 1);
        if ((pt->maxFrstLyrThk) > 0)
          tmpDict1.add("maxFirstLayerThickness", pt->maxFrstLyrThk);
        if ((pt->blNLyr) > 0) tmpDict1.add("nLayers", pt->blNLyr);
        if ((pt->blThkRto) > 0) tmpDict1.add("thicknessRatio", pt->blThkRto);
        patchBoundaryLayers.add(Foam::word("\"" + pt->patchName + "\""),
                                tmpDict1);
      }
      boundaryLayrs.add("patchBoundaryLayers", patchBoundaryLayers);
    }
    meshDict_->add("boundaryLayers", boundaryLayrs);
  }

  // object refinements
  Foam::dictionary objectRefinements("objectRefinements");
  if (!params_->objRefLst.empty()) {
    for (auto ref = (params_->objRefLst).begin();
         ref != (params_->objRefLst).end(); ref++) {
      Foam::dictionary tmpDict1(ref->name);
      for (auto prm = (ref->params).begin(); prm != (ref->params).end();
           prm++) {
        tmpDict1.add(Foam::word(prm->first), Foam::word(prm->second));
      }
      objectRefinements.add(Foam::word(ref->name), tmpDict1);
    }
    meshDict_->add("objectRefinements", objectRefinements);
  }

  // local refinement
  Foam::dictionary localRefinement("localRefinement");
  if (!params_->refPatches.empty()) {
    for (auto pt = (params_->refPatches).begin();
         pt != (params_->refPatches).end(); pt++) {
      Foam::dictionary tmpDict1(pt->patchName);
      if ((pt->aditRefLvls) > 0)
        tmpDict1.add("additionalRefinementLevels", pt->aditRefLvls);
      if ((pt->refThickness) > 0)
        tmpDict1.add("refinementThickness", pt->refThickness);
      if ((pt->cellSize) > 0) tmpDict1.add("cellSize", pt->cellSize);
      localRefinement.add(Foam::word("\"" + pt->patchName + "\""), tmpDict1);
    }
    meshDict_->add("localRefinement", localRefinement);
  }

  // rename boundaries
  Foam::dictionary renameBoundary("renameBoundary");
  if (params_->renBndry.has_value()) {
    const auto &renBndry = params_->renBndry.value();
    renameBoundary.add("defaultName", Foam::word((renBndry).defName));
    renameBoundary.add("defaultType", Foam::word((renBndry).defType));

    Foam::dictionary tmpDict1("newPatchNames");

    for (auto pt = (renBndry).newPatches.begin();
         pt != (renBndry).newPatches.end(); pt++) {
      Foam::dictionary tmpDict2(pt->name);
      tmpDict2.add("newName", Foam::word(pt->newName));
      tmpDict2.add("type", Foam::word(pt->newType));
      tmpDict1.add(Foam::word("\"" + pt->name + "\""), tmpDict2);
    }
    renameBoundary.add("newPatchNames", tmpDict1);
    meshDict_->add("renameBoundary", renameBoundary);
  }

  if (write) {
    // Write meshDict
    Foam::fileName cfmeshDict_ = "system/meshDict";
    Foam::OFstream outcfmeshDict_(cfmeshDict_);
    Foam::IOobject::writeBanner(outcfmeshDict_);
    meshDict_->write(outcfmeshDict_, false);
  }
}
