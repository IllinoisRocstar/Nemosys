#include <cstdlib>

#include <iostream>
#include <string>
#include <vector>

#include "Refinement/AMRFoam.H"
#include "AuxiliaryFunctions.H"

#include <IFstream.H>
#include <IOdictionary.H>
#include <OFstream.H>
#include <ReadFields.H>
#include <Time.H>
#include <argList.H>
#include <cellSet.H>
#include <fvCFD.H>
#include <fvMesh.H>
#include <hexRef8.H>
#include <mapPolyMesh.H>
#include <motionSolver.H>
#include <pointFields.H>
#include <pointMesh.H>
#include <polyMesh.H>
#include <polyTopoChange.H>
#include <sigFpe.H>
#include <surfaceFields.H>
#include <surfaceInterpolate.H>
#include <syncTools.H>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam {

AMRFoam::AMRFoam(const Foam::IOobject& iomesh)
    : dynamicFvMesh(iomesh),
      meshCutter_(*this),
      dumpLevel_(true),
      nRefinementIterations_(0),
      protectedCell_(nCells()) {
  checkForMotion();
  calcProtectedCells();
}

void AMRFoam::checkForMotion() {
  // Check if motion needs to be enabled
  dictionary refineDict(IOdictionary(
      IOobject("dynamicMeshDict", time().constant(), *this,
               IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE, false)));

  word subType = word(refineDict.lookup("dynamicFvMesh"));
  if (subType == "dynamicMotionSolverFvMesh") {
    enableMotion = true;
  }

  if (enableMotion) {
    motionPtr_ = motionSolver::New(*this);
  }
}

void AMRFoam::disableMotion() { enableMotion = false; }

void AMRFoam::calcProtectedCells() {
  const labelList& cellLevel = meshCutter_.cellLevel();
  const labelList& pointLevel = meshCutter_.pointLevel();

  // Set cells that should not be refined.
  // This is currently any cell which does not have 8 anchor points or
  // uses any face which does not have 4 anchor points.
  // Note: do not use cellPoint addressing

  // Count number of points <= cellLevel
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  labelList nAnchors(nCells(), Zero);

  label nProtected = 0;

  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];

    for (const label celli : pCells) {
      if (!protectedCell_.test(celli)) {
        if (pointLevel[pointi] <= cellLevel[celli]) {
          ++nAnchors[celli];

          if (nAnchors[celli] > 8) {
            protectedCell_.set(celli);
            nProtected++;
          }
        }
      }
    }
  }

  // Count number of points <= faceLevel
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Bit tricky since proc face might be one more refined than the owner since
  // the coupled one is refined.

  {
    labelList neiLevel(nFaces());

    for (label facei = 0; facei < nInternalFaces(); ++facei) {
      neiLevel[facei] = cellLevel[faceNeighbour()[facei]];
    }
    for (label facei = nInternalFaces(); facei < nFaces(); ++facei) {
      neiLevel[facei] = cellLevel[faceOwner()[facei]];
    }
    syncTools::swapFaceList(*this, neiLevel);

    bitSet protectedFace(nFaces());

    forAll(faceOwner(), facei) {
      const label faceLevel =
          max(cellLevel[faceOwner()[facei]], neiLevel[facei]);

      const face& f = faces()[facei];

      label nnAnchors = 0;

      for (const label pointi : f) {
        if (pointLevel[pointi] <= faceLevel) {
          ++nnAnchors;

          if (nnAnchors > 4) {
            protectedFace.set(facei);
            break;
          }
        }
      }
    }

    syncTools::syncFaceList(*this, protectedFace, orEqOp<unsigned int>());

    for (label facei = 0; facei < nInternalFaces(); ++facei) {
      if (protectedFace.test(facei)) {
        protectedCell_.set(faceOwner()[facei]);
        nProtected++;
        protectedCell_.set(faceNeighbour()[facei]);
        nProtected++;
      }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); ++facei) {
      if (protectedFace.test(facei)) {
        protectedCell_.set(faceOwner()[facei]);
        nProtected++;
      }
    }

    // Also protect any cells that are less than hex
    forAll(cells(), celli) {
      const cell& cFaces = cells()[celli];

      if (cFaces.size() < 6) {
        if (protectedCell_.set(celli)) {
          nProtected++;
        }

      } else {
        for (const label cfacei : cFaces) {
          if (faces()[cfacei].size() < 4) {
            if (protectedCell_.set(celli)) {
              nProtected++;
            }
            break;
          }
        }
      }
    }

    // Check cells for 8 corner points
    checkEightAnchorPoints(protectedCell_);

    if (enableMotion) {
      // Block boundary cells from refinement
      const polyBoundaryMesh& patches = boundaryMesh();
      wordList ptchTypes = patches.types();

      labelList assignedScore(nCells(), 0);

      forAll(patches, patchI) {
        if (ptchTypes[patchI] != "empty") {
          auto allCells = patches[patchI].faceCells();
          forAll(allCells, cellI) {
            assignedScore[allCells[cellI]] += 1;
            // Find all neighbours and assign the score of one
            auto neiCells = cellCells()[allCells[cellI]];
            forAll(neiCells, nCellI) { assignedScore[neiCells[nCellI]] += 1; }
          }
        }
      }

      forAll(assignedScore, scoreI) {
        if (assignedScore[scoreI] >= 2) {
          protectedCell_.set(scoreI);
        }
      }
    }
  }

  if (!returnReduce(protectedCell_.any(), orOp<bool>())) {
    protectedCell_.clear();
  } else {
    cellSet protectedCells(*this, "protectedCells",
                           HashSetOps::used(protectedCell_));
    Info << "Detected " << returnReduce(protectedCells.size(), sumOp<label>())
         << " cells that are protected from refinement."
         << " Writing these to cellSet " << protectedCells.name() << "."
         << endl;

    protectedCells.write();
  }
}

bool AMRFoam::update() {
  // Re-read dictionary. Chosen since usually -small so trivial amount
  // of time compared to actual refinement. Also very useful to be able
  // to modify on-the-fly.
  dictionary refineDict(
      IOdictionary(IOobject("dynamicMeshDict", time().constant(), *this,
                            IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE,
                            false))
          .optionalSubDict("dynamicRefineMotionFvMeshCoeffs"));

  auto fluxVelocities = refineDict.get<List<Pair<word>>>("correctFluxes");

  // Rework into hashtable.
  correctFluxes_.resize(fluxVelocities.size());
  for (const auto& pr : fluxVelocities) {
    correctFluxes_.insert(pr.first(), pr.second());
  }

  refineDict.readEntry("dumpLevel", dumpLevel_);
  const label refineInterval = refineDict.get<label>("refineInterval");

  bool hasChanged = false;

  if (refineInterval == 0) {
    polyMesh::topoChanging(hasChanged);
    return hasChanged;
  } else if (refineInterval < 0) {
    FatalErrorInFunction << "Illegal refineInterval " << refineInterval << nl
                         << "The refineInterval value should"
                         << " be >= 1." << nl << exit(FatalError);
  }

  if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0) {
    label maxCells = refineDict.get<label>("maxCells");

    if (maxCells < 0) maxCells = 500000;

    if (maxCells <= 0) {
      FatalErrorInFunction
          << "Illegal maximum number of cells " << maxCells << nl
          << "The maxCells setting in the dynamicMeshDict should"
          << " be > 0." << nl << exit(FatalError);
    }

    const label maxRefinement = refineDict.get<label>("maxRefinement");

    if (maxRefinement <= 0) {
      FatalErrorInFunction
          << "Illegal maximum refinement level " << maxRefinement << nl
          << "The maxCells setting in the dynamicMeshDict should"
          << " be > 0." << nl << exit(FatalError);
    }

    const word fieldName(refineDict.get<word>("field"));

    const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

    const scalar lowerRefineLevel = refineDict.get<scalar>("lowerRefineLevel");
    const scalar upperRefineLevel = refineDict.get<scalar>("upperRefineLevel");
    const scalar unrefineBelow =
        refineDict.getOrDefault<scalar>("unrefineLevel", GREAT);
    const scalar unrefineAbove =
        refineDict.getOrDefault<scalar>("unrefineAbove", 50000);
    const label nBufferLayers = refineDict.get<label>("nBufferLayers");

    // Cells marked for refinement or otherwise protected from unrefinement.
    bitSet refineCell(nCells());

    // Determine candidates for refinement (looking at field only)
    selectRefineCandidates(lowerRefineLevel, upperRefineLevel, vFld,
                           refineCell);

    if (globalData().nTotalCells() < maxCells) {
      // Select subset of candidates. Take into account max allowable
      // cells, refinement level, protected cells.
      labelList cellsToRefine(
          selectRefineCells(maxCells, maxRefinement, refineCell));

      const label nCellsToRefine =
          returnReduce(cellsToRefine.size(), sumOp<label>());

      if (nCellsToRefine > 0) {
        // Refine/update mesh and map fields
        autoPtr<mapPolyMesh> map = refine(cellsToRefine);

        // Update refineCell. Note that some of the marked ones have
        // not been refined due to constraints.
        {
          const labelList& cellMap = map().cellMap();
          const labelList& reverseCellMap = map().reverseCellMap();

          bitSet newRefineCell(cellMap.size());

          forAll(cellMap, celli) {
            const label oldCelli = cellMap[celli];
            if ((oldCelli < 0) || (reverseCellMap[oldCelli] != celli) ||
                (refineCell.test(oldCelli))) {
              newRefineCell.set(celli);
            }
          }
          refineCell.transfer(newRefineCell);
        }

        // Extend with a buffer layer to prevent neighbouring points
        // being unrefined.
        for (label i = 0; i < nBufferLayers; ++i) {
          extendMarkedCells(refineCell);
        }

        hasChanged = true;
      }
    }

    {
      // Select unrefineable points that are not marked in refineCell
      labelList pointsToUnrefine(selectUnrefinePoints(
          unrefineAbove, unrefineBelow, refineCell, maxCellField(vFld)));

      const label nSplitPoints =
          returnReduce(pointsToUnrefine.size(), sumOp<label>());

      if (nSplitPoints > 0) {
        // Refine/update mesh
        unrefine(pointsToUnrefine);

        hasChanged = true;
      }
    }

    if ((nRefinementIterations_ % 1) == 0) {
      // Compact refinement history occasionally (how often?).
      // Unrefinement causes holes in the refinementHistory.
      const_cast<refinementHistory&>(meshCutter().history()).compact();
    }
    nRefinementIterations_++;
  }

  writeData();

  writeMesh();

  topoChanging(hasChanged);

  if (hasChanged) {
    // Reset moving flag (if any). If not using inflation we'll not move,
    // if are using inflation any follow on movePoints will set it.
    moving(false);
  } else {
    if (enableMotion) solveMotion();
  }
  return hasChanged;
}

bool AMRFoam::updateAMR(const int& refineInterval, const int& maxRefinement,
                        volScalarField& vFld, const double& lowerRefineLevel,
                        const double& upperRefineLevel,
                        const double& unrefineAbove,
                        const double& unrefineBelow, const int& nBufferLayers,
                        const int& maxCells) {
  bool hasChanged = false;

  if (refineInterval == 0) {
    polyMesh::topoChanging(hasChanged);
    return hasChanged;
  } else if (refineInterval < 0) {
    FatalErrorInFunction << "Illegal refineInterval " << refineInterval << nl
                         << "The refineInterval value should"
                         << " be >= 1." << nl << exit(FatalError);
  }

  if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0) {
    if (maxCells <= 0) {
      FatalErrorInFunction
          << "Illegal maximum number of cells " << maxCells << nl
          << "The maxCells setting in the dynamicMeshDict should"
          << " be > 0." << nl << exit(FatalError);
    }

    if (maxRefinement <= 0) {
      FatalErrorInFunction
          << "Illegal maximum refinement level " << maxRefinement << nl
          << "The maxCells setting in the dynamicMeshDict should"
          << " be > 0." << nl << exit(FatalError);
    }

    // Cells marked for refinement or otherwise protected from unrefinement.
    bitSet refineCell(nCells());

    // Determine candidates for refinement (looking at field only)
    selectRefineCandidates(lowerRefineLevel, upperRefineLevel, vFld,
                           refineCell);

    if (globalData().nTotalCells() < maxCells) {
      // Select subset of candidates. Take into account max allowable
      // cells, refinement level, protected cells.
      labelList cellsToRefine(
          selectRefineCells(maxCells, maxRefinement, refineCell));

      const label nCellsToRefine =
          returnReduce(cellsToRefine.size(), sumOp<label>());

      if (nCellsToRefine > 0) {
        // Refine/update mesh and map fields
        autoPtr<mapPolyMesh> map = refine(cellsToRefine);

        // Update refineCell. Note that some of the marked ones have
        // not been refined due to constraints.
        {
          const labelList& cellMap = map().cellMap();
          const labelList& reverseCellMap = map().reverseCellMap();

          bitSet newRefineCell(cellMap.size());

          forAll(cellMap, celli) {
            const label oldCelli = cellMap[celli];
            if ((oldCelli < 0) || (reverseCellMap[oldCelli] != celli) ||
                (refineCell.test(oldCelli))) {
              newRefineCell.set(celli);
            }
          }
          refineCell.transfer(newRefineCell);
        }

        // Extend with a buffer layer to prevent neighbouring points
        // being unrefined.
        for (label i = 0; i < nBufferLayers; ++i) {
          extendMarkedCells(refineCell);
        }

        hasChanged = true;
      }
    }

    {
      // Select unrefineable points that are not marked in refineCell
      labelList pointsToUnrefine(selectUnrefinePoints(
          unrefineAbove, unrefineBelow, refineCell, maxCellField(vFld)));

      const label nSplitPoints =
          returnReduce(pointsToUnrefine.size(), sumOp<label>());

      if (nSplitPoints > 0) {
        // Refine/update mesh
        unrefine(pointsToUnrefine);
        hasChanged = true;
      }
    }

    if ((nRefinementIterations_ % 10) == 0) {
      // Compact refinement history occasionally (how often?).
      // Unrefinement causes holes in the refinementHistory.
      const_cast<refinementHistory&>(meshCutter().history()).compact();
    }
    nRefinementIterations_++;
  }

  topoChanging(hasChanged);

  if (hasChanged) {
    // Reset moving flag (if any). If not using inflation we'll not move,
    // if are using inflation any follow on movePoints will set it.
    moving(false);
  } else {
    if (enableMotion) solveMotion();
  }

  if (writeField) vFld.write();

  if (writeRefHistory) writeData();

  if (writeMeshData) writeMesh();

  return hasChanged;
}

// New method for machine learning model
bool AMRFoam::updateAMRML(const int& refineInterval, const int& maxRefinement,
                          const int& nBufferLayers, const int& maxCells,
                          volScalarField& vFld) {
  bool hasChanged = false;
  if (refineInterval == 0) {
    polyMesh::topoChanging(hasChanged);
    return hasChanged;
  } else if (refineInterval < 0) {
    FatalErrorInFunction << "Illegal refineInterval " << refineInterval << nl
                         << "The refineInterval value should"
                         << " be >= 1." << nl << exit(FatalError);
  }

  // if (maxCells < 0)
  //   maxCells = 500000;

  // Cannot refine at time = 0;
  if (fvMesh::time().timeIndex() > 0 &&
      fvMesh::time().timeIndex() % refineInterval == 0) {
    if (maxCells <= 0) {
      FatalErrorInFunction << "Illegal maximum number of cells " << maxCells
                           << nl << "The maxCells value should"
                           << " be > 0." << nl << exit(FatalError);
    }

    if (maxRefinement <= 0) {
      FatalErrorInFunction << "Illegal maximum refinement level "
                           << maxRefinement << nl << "The maxCells value should"
                           << " be > 0." << nl << exit(FatalError);
    }

    // Cells marked for refinement or otherwise protected from unrefinement.
    bitSet refineCell(primitiveMesh::nCells());

    // Selecting refinement candidates based on lower and upper refinement
    // levels in scalar field
    selectRefineCandidates(0.5, 1.5, vFld, refineCell);

    // Refinement Procedure Loop
    if (polyMesh::globalData().nTotalCells() < maxCells) {
      // Select subset of candidates. Take into account max allowable
      // cells, refinement level, protected cells.
      labelList cellsToRefine(
          selectRefineCells(maxCells, maxRefinement, refineCell));

      label nCellsToRefine = returnReduce(cellsToRefine.size(), sumOp<label>());

      if (nCellsToRefine > 0) {
        // Refine/update mesh and map fields
        autoPtr<mapPolyMesh> map = refine(cellsToRefine);

        // Update refineCell. Note that some of the marked ones have not been
        // refined due to constraints.
        const labelList& cellMap = map().cellMap();
        const labelList& reverseCellMap = map().reverseCellMap();

        bitSet newRefineCell(cellMap.size());

        forAll(cellMap, celli) {
          label oldCelli = cellMap[celli];

          if (oldCelli < 0)
            newRefineCell.set(celli, 1);
          else if (reverseCellMap[oldCelli] != celli)
            newRefineCell.set(celli, 1);
          else
            newRefineCell.set(celli, refineCell.get(oldCelli));
        }
        refineCell.transfer(newRefineCell);

        // Extended with a buffer layer to prevent neighbouring points
        // being unrefined
        for (label i = 0; i < nBufferLayers; i++) {
          extendMarkedCells(refineCell);
        }

        hasChanged = true;
      }
    }

    // Select unrefinable points that are not marked in refineCell
    labelList pointsToUnrefine(
        selectUnrefinePoints(-1, 0.5, refineCell, maxCellField(vFld)));

    label nSplitPoints = returnReduce(pointsToUnrefine.size(), sumOp<label>());
    autoPtr<mapPolyMesh> mapUnrefine;

    if (nSplitPoints > 0) {
      // Refine/update mesh
      unrefine(pointsToUnrefine);
      hasChanged = true;
    }

    if ((nRefinementIterations_ % 10) == 0) {
      // Unrefinement causes holes in the refinementHistory
      const_cast<refinementHistory&>(meshCutter().history()).compact();
    }
    nRefinementIterations_++;
  }

  polyMesh::topoChanging(hasChanged);
  if (hasChanged) {
    polyMesh::moving(false);
  }

  return hasChanged;
}

labelList AMRFoam::selectUnrefinePoints(const scalar unrefineAbove,
                                        const scalar unrefineBelow,
                                        const bitSet& markedCell,
                                        const scalarField& pFld) const {
  // All points that can be unrefined
  const labelList splitPoints(meshCutter_.getSplitPoints());

  const labelListList& pointCells = this->pointCells();

  // If we have any protected cells make sure they also are not being
  // unrefined

  bitSet protectedPoint(nPoints());

  if (!protectedCell_.empty()) {
    // Get all points on a protected cell
    forAll(pointCells, pointi) {
      for (const label celli : pointCells[pointi]) {
        if (protectedCell_.test(celli)) {
          protectedPoint.set(pointi);
          break;
        }
      }
    }

    syncTools::syncPointList(*this, protectedPoint, orEqOp<unsigned int>(), 0u);

    DebugInfo << "From " << returnReduce(protectedCell_.count(), sumOp<label>())
              << " protected cells found "
              << returnReduce(protectedPoint.count(), sumOp<label>())
              << " protected points." << endl;
  }

  DynamicList<label> newSplitPoints(splitPoints.size());

  for (const label pointi : splitPoints) {
    if ((!protectedPoint[pointi] && pFld[pointi] < unrefineBelow) ||
        (!protectedPoint[pointi] && pFld[pointi] > unrefineAbove)) {
      // Check that all cells are not marked
      bool hasMarked = false;

      for (const label celli : pointCells[pointi]) {
        if (markedCell.test(celli)) {
          hasMarked = true;
          break;
        }
      }

      if (!hasMarked) {
        newSplitPoints.append(pointi);
      }
    }
  }

  newSplitPoints.shrink();

  // Guarantee 2:1 refinement after unrefinement
  labelList consistentSet(
      meshCutter_.consistentUnrefinement(newSplitPoints, false));
  Info << "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
       << " split points out of a possible "
       << returnReduce(splitPoints.size(), sumOp<label>()) << "." << endl;

  return consistentSet;
}

scalarField AMRFoam::cell2Pt(const scalarField& vFld) {
  scalarField pFld(nPoints());

  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];

    scalar sum = 0.0;
    forAll(pCells, i) sum += vFld[pCells[i]];

    pFld[pointi] = sum / pCells.size();
  }

  return pFld;
}

void AMRFoam::writeData() {
  // Just write refiment history in current time
  const_cast<hexRef8&>(meshCutter_).setInstance(time().timeName());
  meshCutter_.write(true);
}

const polyMesh& AMRFoam::getMesh() { return meshCutter_.mesh(); }

void AMRFoam::writeMesh() {
  auto strmOpt =
      IOstreamOption(time().writeFormat(), time().writeCompression());
  dynamicFvMesh::writeObject(strmOpt, true);
}

void AMRFoam::transformField(volScalarField& inField) {
  double maxNum = 0;
  for (int i = 0; i < inField.size(); i++) {
    if (inField[i] < 0) {
      // Do Nothing
    } else if (inField[i] > maxNum) {
      maxNum = inField[i];
    } else {
      // Do Nothing
    }
  }

  for (int i = 0; i < inField.size(); i++) inField[i] = inField[i] / maxNum;
}

volScalarField AMRFoam::readIncomingCellField(const std::string& inName,
                                              const std::string& outName) {
  std::cout << " - Parsing Input file ... " << std::endl;
  std::ifstream inFld(inName);
  if (inFld.is_open()) {
    std::string linesOut;              // Temporary string for multiple uses
    std::vector<std::string> myLines;  // Stores whole file line by line
    while (std::getline(inFld, linesOut)) myLines.push_back(linesOut);

    incomingField.clear();
    incomingField.resize(myLines.size());

    for (int i = 0; i < (int)myLines.size(); i++)
      incomingField[i] =
          std::strtod(nemAux::strToChar(myLines[i]).get(), nullptr);

    volScalarField returnFld = initScalarField(outName);

    for (int i = 0; i < returnFld.size(); i++) returnFld[i] = incomingField[i];

    return returnFld;
  } else {
    std::cerr << "Cannot open/find " << inName << " file!" << std::endl;
    throw;
  }
}

volScalarField AMRFoam::assignToVolScalarField(const std::vector<int>& vec) {
  volScalarField vFld = initScalarField("vFld");
  for (int i = 0; i < (int)vec.size(); i++) vFld[i] = vec[i];
  return vFld;
}

pointScalarField AMRFoam::readIncomingPtField(const std::string& inName,
                                              const std::string& outName) {
  std::cout << " - Parsing Input file ... " << std::endl;
  std::ifstream inFld(inName);
  if (inFld.is_open()) {
    std::string linesOut;              // Temporary string for multiple uses
    std::vector<std::string> myLines;  // Stores whole file line by line
    while (std::getline(inFld, linesOut)) myLines.push_back(linesOut);

    incomingField.clear();
    incomingField.resize(myLines.size());

    for (int i = 0; i < (int)myLines.size(); i++) {
      incomingField[i] =
          std::strtod(nemAux::strToChar(myLines[i]).get(), nullptr);
    }

    pointScalarField returnFld = initPntSField(outName);

    for (int i = 0; i < returnFld.size(); i++) returnFld[i] = incomingField[i];

    return returnFld;
  } else {
    std::cerr << "Cannot open/find " << inName << " file!" << std::endl;
    throw;
  }
}

volScalarField AMRFoam::readInitialField(const std::string& fldName) {
  volScalarField meshFieldXY(IOobject(fldName, time().timeName(), *this,
                                      IOobject::NO_READ, IOobject::AUTO_WRITE),
                             *this);

  return meshFieldXY;
}

volScalarField AMRFoam::pF2vF(const pointScalarField& pntF,
                              const std::string& outName) {
  volScalarField returnFld = initScalarField(outName);
  return returnFld;
}

volScalarField AMRFoam::initScalarField(const std::string& fldName) {
  volScalarField emptyScalar(IOobject(fldName, time().timeName(), *this,
                                      IOobject::NO_READ, IOobject::AUTO_WRITE),
                             *this, dimensionedScalar("scalar", dimLength, 0));

  return emptyScalar;
}

volScalarField AMRFoam::initGradientField(const std::string& fldName) {
  volScalarField emptyScalar(IOobject(fldName, time().timeName(), *this,
                                      IOobject::NO_READ, IOobject::AUTO_WRITE),
                             *this, dimensionedScalar("scalar", dimless, 0));

  return emptyScalar;
}

pointScalarField AMRFoam::initPntSField(const std::string& fldName) {
  // Initialize a pointMesh object
  pointMesh pMesh(*this);

  // Zero pointScalarField declaration. InternalField contains all mesh
  // points. Can be overwritten with incoming stream.
  pointScalarField initFieldZero(
      IOobject(fldName, time().timeName(), *this, IOobject::NO_READ,
               IOobject::AUTO_WRITE),
      pMesh, dimensionedScalar("scalar", dimLength, 0));

  return initFieldZero;
}

void AMRFoam::enableRefHistoryData() { writeRefHistory = true; }

void AMRFoam::enableUpdatedField() { writeField = true; }

void AMRFoam::enableMeshWriting() { writeMeshData = true; }

volScalarField AMRFoam::getGradient(const volScalarField& fldGrd) {
  volVectorField returnGrad(IOobject("fldName", time().timeName(), *this,
                                     IOobject::NO_READ, IOobject::NO_WRITE),
                            *this, dimensionedVector("vector", dimless, Zero));

  returnGrad = fvc::grad(fldGrd);
  returnGrad.write();
  volScalarField returnScalar = initGradientField("gradient");
  returnScalar = mag(returnGrad);

  return returnScalar;
}

void AMRFoam::calculateProtectedCells(bitSet& unrefineableCell) const {
  if (protectedCell_.empty()) {
    unrefineableCell.clear();
    return;
  }

  const labelList& cellLevel = meshCutter_.cellLevel();

  unrefineableCell = protectedCell_;

  // Get neighbouring cell level
  labelList neiLevel(nFaces() - nInternalFaces());

  for (label facei = nInternalFaces(); facei < nFaces(); ++facei) {
    neiLevel[facei - nInternalFaces()] = cellLevel[faceOwner()[facei]];
  }
  syncTools::swapBoundaryFaceList(*this, neiLevel);

  bitSet seedFace;

  while (true) {
    // Pick up faces on border of protected cells
    seedFace.reset();
    seedFace.resize(nFaces());

    for (label facei = 0; facei < nInternalFaces(); ++facei) {
      const label own = faceOwner()[facei];
      const label nei = faceNeighbour()[facei];

      if (
          // Protected owner
          (unrefineableCell.test(own) && (cellLevel[nei] > cellLevel[own])) ||
          // Protected neighbour
          (unrefineableCell.test(nei) && (cellLevel[own] > cellLevel[nei]))) {
        seedFace.set(facei);
      }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); facei++) {
      const label own = faceOwner()[facei];

      if (
          // Protected owner
          (unrefineableCell.test(own) &&
           (neiLevel[facei - nInternalFaces()] > cellLevel[own]))) {
        seedFace.set(facei);
      }
    }

    syncTools::syncFaceList(*this, seedFace, orEqOp<unsigned int>());

    // Extend unrefineableCell
    bool hasExtended = false;

    for (label facei = 0; facei < nInternalFaces(); ++facei) {
      if (seedFace.test(facei)) {
        if (unrefineableCell.set(faceOwner()[facei])) {
          hasExtended = true;
        }
        if (unrefineableCell.set(faceNeighbour()[facei])) {
          hasExtended = true;
        }
      }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); ++facei) {
      if (seedFace.test(facei)) {
        const label own = faceOwner()[facei];

        if (unrefineableCell.set(own)) {
          hasExtended = true;
        }
      }
    }

    if (!returnReduce(hasExtended, orOp<bool>())) {
      break;
    }
  }
}

void AMRFoam::readDict() {
  dictionary refineDict(
      IOdictionary(IOobject("dynamicMeshDict", time().constant(), *this,
                            IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE,
                            false))
          .optionalSubDict("dynamicRefineMotionFvMeshCoeffs"));

  List<Pair<word>> fluxVelocities =
      List<Pair<word>>(refineDict.lookup("correctFluxes"));
  // Rework into hashtable.
  correctFluxes_.resize(fluxVelocities.size());
  forAll(fluxVelocities, i)
      correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);

  dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
}

void AMRFoam::mapFields(const mapPolyMesh& mpm) {
  dynamicFvMesh::mapFields(mpm);

  // Correct the flux for modified/added faces. All the faces which only
  // have been renumbered will already have been handled by the mapping.
  {
    const labelList& faceMap = mpm.faceMap();
    const labelList& reverseFaceMap = mpm.reverseFaceMap();

    // Storage for any master faces. These will be the original faces
    // on the coarse cell that get split into four (or rather the
    // master face gets modified and three faces get added from the master)
    // Estimate number of faces created

    bitSet masterFaces(nFaces());

    forAll(faceMap, facei) {
      const label oldFacei = faceMap[facei];

      if (oldFacei >= 0) {
        const label masterFacei = reverseFaceMap[oldFacei];

        if (masterFacei < 0) {
          FatalErrorInFunction << "Problem: should not have removed faces"
                               << " when refining." << nl << "face:" << facei
                               << endl
                               << abort(FatalError);
        } else if (masterFacei != facei) {
          masterFaces.set(masterFacei);
        }
      }
    }

    if (debug) {
      Pout << "Found " << masterFaces.count() << " split faces " << endl;
    }

    HashTable<surfaceScalarField*> fluxes(lookupClass<surfaceScalarField>());
    forAllIters(fluxes, iter) {
      if (!correctFluxes_.found(iter.key())) {
        WarningInFunction
            << "Cannot find surfaceScalarField " << iter.key()
            << " in user-provided flux mapping table " << correctFluxes_ << endl
            << "    The flux mapping table is used to recreate the"
            << " flux on newly created faces." << endl
            << "    Either add the entry if it is a flux or use (" << iter.key()
            << " none) to suppress this warning." << endl;
        continue;
      }

      const word& UName = correctFluxes_[iter.key()];

      if (UName == "none") {
        continue;
      }

      surfaceScalarField& phi = *iter();

      if (UName == "NaN") {
        Pout << "Setting surfaceScalarField " << iter.key() << " to NaN"
             << endl;

        sigFpe::fillNan(phi.primitiveFieldRef());

        continue;
      }

      if (debug) {
        Pout << "Mapping flux " << iter.key() << " using interpolated flux "
             << UName << endl;
      }

      const surfaceScalarField phiU(
          fvc::interpolate(lookupObject<volVectorField>(UName)) & Sf());

      // Recalculate new internal faces.
      for (label facei = 0; facei < nInternalFaces(); ++facei) {
        const label oldFacei = faceMap[facei];

        if (oldFacei == -1) {
          // Inflated/appended
          phi[facei] = phiU[facei];
        } else if (reverseFaceMap[oldFacei] != facei) {
          // face-from-masterface
          phi[facei] = phiU[facei];
        }
      }

      // Recalculate new boundary faces.
      surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

      forAll(phiBf, patchi) {
        fvsPatchScalarField& patchPhi = phiBf[patchi];
        const fvsPatchScalarField& patchPhiU = phiU.boundaryField()[patchi];

        label facei = patchPhi.patch().start();

        forAll(patchPhi, i) {
          const label oldFacei = faceMap[facei];

          if (oldFacei == -1) {
            // Inflated/appended
            patchPhi[i] = patchPhiU[i];
          } else if (reverseFaceMap[oldFacei] != facei) {
            // face-from-masterface
            patchPhi[i] = patchPhiU[i];
          }

          ++facei;
        }
      }

      // Update master faces
      for (const label facei : masterFaces) {
        if (isInternalFace(facei)) {
          phi[facei] = phiU[facei];
        } else {
          const label patchi = boundaryMesh().whichPatch(facei);
          const label i = facei - boundaryMesh()[patchi].start();

          const fvsPatchScalarField& patchPhiU = phiU.boundaryField()[patchi];

          fvsPatchScalarField& patchPhi = phiBf[patchi];

          patchPhi[i] = patchPhiU[i];
        }
      }
    }
  }

  // Correct the flux for injected faces - these are the faces which have
  // no correspondence to the old mesh (i.e. added without a masterFace, edge
  // or point). An example is the internal faces from hexRef8.
  {
    const labelList& faceMap = mpm.faceMap();

    mapNewInternalFaces<scalar>(this->Sf(), this->magSf(), faceMap);
    mapNewInternalFaces<vector>(this->Sf(), this->magSf(), faceMap);

    // No oriented fields of more complex type
    mapNewInternalFaces<sphericalTensor>(faceMap);
    mapNewInternalFaces<symmTensor>(faceMap);
    mapNewInternalFaces<tensor>(faceMap);
  }
}

// Refines cells, maps fields and recalculates (an approximate) flux
autoPtr<mapPolyMesh> Foam::AMRFoam::refine(const labelList& cellsToRefine) {
  // Mesh changing engine.
  polyTopoChange meshMod(*this);

  // Play refinement commands into mesh changer.
  meshCutter_.setRefinement(cellsToRefine, meshMod);

  // Create mesh (with inflation), return map from old to new mesh.
  // autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
  autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

  Info << "Refined from " << returnReduce(map().nOldCells(), sumOp<label>())
       << " to " << globalData().nTotalCells() << " cells." << endl;

  if (debug) {
    // Check map.
    for (label facei = 0; facei < nInternalFaces(); ++facei) {
      const label oldFacei = map().faceMap()[facei];

      if (oldFacei >= nInternalFaces()) {
        FatalErrorInFunction << "New internal face:" << facei
                             << " fc:" << faceCentres()[facei]
                             << " originates from boundary oldFace:" << oldFacei
                             << abort(FatalError);
      }
    }
  }

  // Update fields
  updateMesh(*map);

  if (enableMotion) {
    // Solve for mesh motion and move points
    motionPtr_->updateMesh(*map);
    movePoints(motionPtr_->newPoints());
  }

  // Update numbering of cells/vertices.
  meshCutter_.updateMesh(*map);

  // Update numbering of protectedCell_
  if (!protectedCell_.empty()) {
    bitSet newProtectedCell(nCells());

    forAll(newProtectedCell, celli) {
      const label oldCelli = map().cellMap()[celli];
      if (protectedCell_.test(oldCelli)) {
        newProtectedCell.set(celli);
      }
    }
    protectedCell_.transfer(newProtectedCell);
  }

  // Debug: Check refinement levels (across faces only)
  meshCutter_.checkRefinementLevels(-1, labelList());

  return map;
}

autoPtr<mapPolyMesh> Foam::AMRFoam::unrefine(const labelList& splitPoints) {
  polyTopoChange meshMod(*this);

  // Play refinement commands into mesh changer.
  meshCutter_.setUnrefinement(splitPoints, meshMod);

  // Save information on faces that will be combined
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // Find the faceMidPoints on cells to be combined.
  // for each face resulting of split of face into four store the
  // midpoint
  Map<label> faceToSplitPoint(3 * splitPoints.size());

  {
    for (const label pointi : splitPoints) {
      const labelList& pEdges = pointEdges()[pointi];

      for (const label edgei : pEdges) {
        const label otherPointi = edges()[edgei].otherVertex(pointi);

        const labelList& pFaces = pointFaces()[otherPointi];

        for (const label facei : pFaces) {
          faceToSplitPoint.insert(facei, otherPointi);
        }
      }
    }
  }

  // Change mesh and generate map.
  // autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
  autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

  Info << "Unrefined from " << returnReduce(map().nOldCells(), sumOp<label>())
       << " to " << globalData().nTotalCells() << " cells." << endl;

  // Update fields
  updateMesh(*map);

  if (enableMotion) {
    // Solve for mesh motion and move points
    motionPtr_->updateMesh(*map);
    movePoints(motionPtr_->newPoints());
  }

  // Correct the flux for modified faces.
  {
    const labelList& reversePointMap = map().reversePointMap();
    const labelList& reverseFaceMap = map().reverseFaceMap();

    HashTable<surfaceScalarField*> fluxes(lookupClass<surfaceScalarField>());
    forAllIters(fluxes, iter) {
      if (!correctFluxes_.found(iter.key())) {
        WarningInFunction
            << "Cannot find surfaceScalarField " << iter.key()
            << " in user-provided flux mapping table " << correctFluxes_ << endl
            << "    The flux mapping table is used to recreate the"
            << " flux on newly created faces." << endl
            << "    Either add the entry if it is a flux or use (" << iter.key()
            << " none) to suppress this warning." << endl;
        continue;
      }

      const word& UName = correctFluxes_[iter.key()];

      if (UName == "none") {
        continue;
      }

      DebugInfo << "Mapping flux " << iter.key() << " using interpolated flux "
                << UName << endl;

      surfaceScalarField& phi = *iter();
      surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

      const surfaceScalarField phiU(
          fvc::interpolate(lookupObject<volVectorField>(UName)) & Sf());

      forAllConstIters(faceToSplitPoint, iter2) {
        const label oldFacei = iter2.key();
        const label oldPointi = iter2.val();
        if (reversePointMap[oldPointi] < 0) {
          // midpoint was removed. See if face still exists.
          const label facei = reverseFaceMap[oldFacei];

          if (facei >= 0) {
            if (isInternalFace(facei)) {
              phi[facei] = phiU[facei];
            } else {
              label patchi = boundaryMesh().whichPatch(facei);
              label i = facei - boundaryMesh()[patchi].start();

              const fvsPatchScalarField& patchPhiU =
                  phiU.boundaryField()[patchi];
              fvsPatchScalarField& patchPhi = phiBf[patchi];
              patchPhi[i] = patchPhiU[i];
            }
          }
        }
      }
    }
  }

  // Update numbering of cells/vertices.
  meshCutter_.updateMesh(*map);

  // Update numbering of protectedCell_
  if (!protectedCell_.empty()) {
    bitSet newProtectedCell(nCells());

    forAll(newProtectedCell, celli) {
      const label oldCelli = map().cellMap()[celli];
      if (protectedCell_.test(oldCelli)) {
        newProtectedCell.set(celli);
      }
    }
    protectedCell_.transfer(newProtectedCell);
  }

  // Debug: Check refinement levels (across faces only)
  meshCutter_.checkRefinementLevels(-1, labelList());

  return map;
}

scalarField AMRFoam::maxPointField(const scalarField& pFld) const {
  scalarField vFld(nCells(), -GREAT);

  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];

    for (const label celli : pCells) {
      vFld[celli] = max(vFld[celli], pFld[pointi]);
    }
  }
  return vFld;
}

scalarField AMRFoam::maxCellField(const volScalarField& vFld) const {
  scalarField pFld(nPoints(), -GREAT);

  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];
    for (const label celli : pCells) {
      pFld[pointi] = max(pFld[pointi], vFld[celli]);
    }
  }
  return pFld;
}

scalarField AMRFoam::cellToPoint(const scalarField& vFld) const {
  scalarField pFld(nPoints());

  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];

    scalar sum = 0.0;
    for (const label celli : pCells) {
      sum += vFld[celli];
    }
    pFld[pointi] = sum / pCells.size();
  }
  return pFld;
}

scalarField AMRFoam::error(const scalarField& fld, const scalar minLevel,
                           const scalar maxLevel) const {
  scalarField c(fld.size(), scalar(-1));

  forAll(fld, i) {
    scalar err = min(fld[i] - minLevel, maxLevel - fld[i]);

    if (err >= 0) {
      c[i] = err;
    }
  }
  return c;
}

void AMRFoam::selectRefineCandidates(const scalar lowerRefineLevel,
                                     const scalar upperRefineLevel,
                                     const scalarField& vFld,
                                     bitSet& candidateCell) const {
  // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
  // higher more desirable to be refined).
  scalarField cellError(maxPointField(
      error(cellToPoint(vFld), lowerRefineLevel, upperRefineLevel)));

  // Mark cells that are candidates for refinement.
  forAll(cellError, celli) {
    if (cellError[celli] > 0) {
      candidateCell.set(celli);
    }
  }
}

labelList AMRFoam::selectRefineCells(const label maxCells,
                                     const label maxRefinement,
                                     const bitSet& candidateCell) const {
  // Every refined cell causes 7 extra cells
  label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

  const labelList& cellLevel = meshCutter_.cellLevel();

  // Mark cells that cannot be refined since they would trigger refinement
  // of protected cells (since 2:1 cascade)
  bitSet unrefineableCell;
  calculateProtectedCells(unrefineableCell);

  // Count current selection
  auto nLocalCandidates = (label)candidateCell.count();
  label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

  // Collect all cells
  DynamicList<label> candidates(nLocalCandidates);

  if (nCandidates < nTotToRefine) {
    for (const label celli : candidateCell) {
      if ((!unrefineableCell.test(celli)) && cellLevel[celli] < maxRefinement) {
        candidates.append(celli);
      }
    }
  } else {
    // Sort by error? For now just truncate.
    for (label level = 0; level < maxRefinement; ++level) {
      for (const label celli : candidateCell) {
        if ((!unrefineableCell.test(celli)) && cellLevel[celli] == level) {
          candidates.append(celli);
        }
      }

      if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine) {
        break;
      }
    }
  }

  // Guarantee 2:1 refinement after refinement
  labelList consistentSet(meshCutter_.consistentRefinement(
      candidates.shrink(),
      true  // Add to set to guarantee 2:1
      ));

  Info << "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
       << " cells for refinement out of " << globalData().nTotalCells() << "."
       << endl;

  return consistentSet;
}

void AMRFoam::extendMarkedCells(bitSet& markedCell) const {
  // Mark faces using any marked cell
  bitSet markedFace(nFaces());

  for (const label celli : markedCell) {
    markedFace.set(cells()[celli]);  // set multiple faces
  }

  syncTools::syncFaceList(*this, markedFace, orEqOp<unsigned int>());

  // Update cells using any markedFace
  for (label facei = 0; facei < nInternalFaces(); ++facei) {
    if (markedFace.test(facei)) {
      markedCell.set(faceOwner()[facei]);
      markedCell.set(faceNeighbour()[facei]);
    }
  }
  for (label facei = nInternalFaces(); facei < nFaces(); ++facei) {
    if (markedFace.test(facei)) {
      markedCell.set(faceOwner()[facei]);
    }
  }
}

void AMRFoam::checkEightAnchorPoints(bitSet& protectedCell) const {
  const labelList& cellLevel = meshCutter_.cellLevel();
  const labelList& pointLevel = meshCutter_.pointLevel();

  labelList nAnchorPoints(nCells(), Zero);

  forAll(pointLevel, pointi) {
    const labelList& pCells = pointCells(pointi);

    for (const label celli : pCells) {
      if (pointLevel[pointi] <= cellLevel[celli]) {
        // Check if cell has already 8 anchor points -> protect cell
        if (nAnchorPoints[celli] == 8) {
          protectedCell.set(celli);
        }

        if (!protectedCell.test(celli)) {
          ++nAnchorPoints[celli];
        }
      }
    }
  }

  forAll(protectedCell, celli) {
    if (nAnchorPoints[celli] != 8) {
      protectedCell.set(celli);
    }
  }
}

bool AMRFoam::writeObject(IOstream::streamFormat fmt,
                          IOstream::versionNumber ver,
                          IOstream::compressionType cmp,
                          const bool valid) const {
  // Force refinement data to go to the current time directory.
  const_cast<hexRef8&>(meshCutter_).setInstance(time().timeName());
  auto strmOpt =
      IOstreamOption(time().writeFormat(), time().writeCompression());
  bool writeOk =
      (dynamicFvMesh::writeObject(strmOpt, true) && meshCutter_.write(valid));

  if (dumpLevel_) {
    volScalarField scalarCellLevel(
        IOobject("cellLevel", time().timeName(), *this, IOobject::NO_READ,
                 IOobject::AUTO_WRITE, false),
        *this, dimensionedScalar("level", dimless, 0));

    const labelList& cellLevel = meshCutter_.cellLevel();

    forAll(cellLevel, celli) { scalarCellLevel[celli] = cellLevel[celli]; }

    writeOk = writeOk && scalarCellLevel.write();
  }

  return writeOk;
}

void AMRFoam::solveMotion() {
  movePoints(motionPtr_->newPoints());

  if (foundObject<volVectorField>("U")) {
    lookupObjectRef<volVectorField>("U").correctBoundaryConditions();
  }
}

template <class T>
void AMRFoam::mapNewInternalFaces(
    const labelList& faceMap,
    GeometricField<T, fvsPatchField, surfaceMesh>& sFld) {
  typedef GeometricField<T, fvsPatchField, surfaceMesh> GeoField;

  //- Make flat field for ease of looping
  Field<T> tsFld(this->nFaces(), Zero);
  SubField<T>(tsFld, this->nInternalFaces()) = sFld.internalField();

  const typename GeoField::Boundary& bFld = sFld.boundaryField();
  forAll(bFld, patchi) {
    label facei = this->boundaryMesh()[patchi].start();
    for (const T& val : bFld[patchi]) {
      tsFld[facei++] = val;
    }
  }

  const labelUList& owner = this->faceOwner();
  const labelUList& neighbour = this->faceNeighbour();
  const cellList& cells = this->cells();

  for (label facei = 0; facei < nInternalFaces(); facei++) {
    label oldFacei = faceMap[facei];

    // Map surface field on newly generated faces by obtaining the
    // hull of the outside faces
    if (oldFacei == -1) {
      // Loop over all owner/neighbour cell faces
      // and find already mapped ones (master-faces):
      T tmpValue(pTraits<T>::zero);
      label counter = 0;

      const cell& ownFaces = cells[owner[facei]];
      for (auto ownFacei : ownFaces) {
        if (faceMap[ownFacei] != -1) {
          tmpValue += tsFld[ownFacei];
          counter++;
        }
      }

      const cell& neiFaces = cells[neighbour[facei]];
      for (auto neiFacei : neiFaces) {
        if (faceMap[neiFacei] != -1) {
          tmpValue += tsFld[neiFacei];
          counter++;
        }
      }

      if (counter > 0) {
        sFld[facei] = tmpValue / counter;
      }
    }
  }
}

template <class T>
void AMRFoam::mapNewInternalFaces(const labelList& faceMap) {
  typedef GeometricField<T, fvsPatchField, surfaceMesh> GeoField;
  HashTable<GeoField*> sFlds(this->objectRegistry::lookupClass<GeoField>());

  forAllIters(sFlds, iter) {
    // if (mapSurfaceFields_.found(iter.key()))
    {
      DebugInfo << "dynamicRefineMotionFvMesh::mapNewInternalFaces():"
                << " Mapping new internal faces by interpolation on "
                << iter.key() << endl;

      GeoField& sFld = *iter();

      if (sFld.oriented()()) {
        WarningInFunction << "Ignoring mapping oriented field " << sFld.name()
                          << " since of type " << sFld.type() << endl;
      } else {
        mapNewInternalFaces(faceMap, sFld);
      }
    }
  }
}

template <class T>
void AMRFoam::mapNewInternalFaces(const surfaceVectorField& Sf,
                                  const surfaceScalarField& magSf,
                                  const labelList& faceMap) {
  typedef GeometricField<T, fvsPatchField, surfaceMesh> GeoField;
  HashTable<GeoField*> sFlds(this->objectRegistry::lookupClass<GeoField>());

  forAllIters(sFlds, iter) {
    // if (mapSurfaceFields_.found(iter.key()))
    {
      DebugInfo << "dynamicRefineMotionFvMesh::mapNewInternalFaces():"
                << " Mapping new internal faces by interpolation on "
                << iter.key() << endl;

      GeoField& sFld = *iter();

      if (sFld.oriented()()) {
        DebugInfo << "dynamicRefineMotionFvMesh::mapNewInternalFaces(): "
                  << "Converting oriented field " << iter.key()
                  << " to intensive field and mapping" << endl;

        // Assume any oriented field is face area weighted (i.e. a flux)
        // Convert to intensive (& oriented) before mapping. Untested.

        typedef GeometricField<typename outerProduct<vector, T>::type,
                               fvsPatchField, surfaceMesh>
            NormalGeoField;

        // Convert to intensive and non oriented
        NormalGeoField fFld(sFld * Sf / Foam::sqr(magSf));

        // Interpolate
        mapNewInternalFaces(faceMap, fFld);

        // Convert back to extensive and oriented
        sFld = (fFld & Sf);
      } else {
        mapNewInternalFaces(faceMap, sFld);
      }
    }
  }
}

}  // namespace Foam

// ************************************************************************* //