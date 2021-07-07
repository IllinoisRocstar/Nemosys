#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "Refinement/AMRFoam.H"
#include "AuxiliaryFunctions.H"

#include <IFstream.H>
#include <IOdictionary.H>
#include <OFstream.H>
#include <ReadFields.H>
#include <Time.H>
#include <addToRunTimeSelectionTable.H>
#include <argList.H>
#include <cellSet.H>
#include <fvCFD.H>
#include <fvMesh.H>
#include <hexRef8.H>
#include <labelIOList.H>
#include <mapPolyMesh.H>
#include <motionSolver.H>
#include <multiDirRefinement.H>
#include <pointFields.H>
#include <pointMesh.H>
#include <pointPatchField.H>
#include <polyMesh.H>
#include <polyTopoChange.H>
#include <sigFpe.H>
#include <surfaceFields.H>
#include <surfaceInterpolate.H>
#include <syncTools.H>
#include <timeSelector.H>
#include <vectorIOField.H>
#include <volFields.H>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam {

AMRFoam::AMRFoam(const Foam::IOobject& iomesh)
    : dynamicFvMesh(iomesh),
      // motionPtr_(motionSolver::New(*this)),
      meshCutter_(*this),
      dumpLevel_(false),
      nRefinementIterations_(0),
      protectedCell_(nCells(), 0) {
  // Nothing
  const labelList& cellLevel = meshCutter_.cellLevel();
  const labelList& pointLevel = meshCutter_.pointLevel();

  // Set cells that should not be refined.
  // This is currently any cell which does not have 8 anchor points or
  // uses any face which does not have 4 anchor points.
  // Note: do not use cellPoint addressing

  // Count number of points <= cellLevel
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  labelList nAnchors(nCells(), 0);

  label nProtected = 0;

  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];
    forAll(pCells, i) {
      label celli = pCells[i];
      if (!protectedCell_.get(celli)) {
        if (pointLevel[pointi] <= cellLevel[celli]) {
          nAnchors[celli]++;
          if (nAnchors[celli] > 8) {
            protectedCell_.set(celli, 1);
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

    for (label facei = 0; facei < nInternalFaces(); facei++)
      neiLevel[facei] = cellLevel[faceNeighbour()[facei]];

    for (label facei = nInternalFaces(); facei < nFaces(); facei++)
      neiLevel[facei] = cellLevel[faceOwner()[facei]];

    syncTools::swapFaceList(*this, neiLevel);

    boolList protectedFace(nFaces(), false);

    forAll(faceOwner(), facei) {
      label faceLevel = max(cellLevel[faceOwner()[facei]], neiLevel[facei]);

      const face& f = faces()[facei];

      label nAnchors = 0;

      forAll(f, fp) {
        if (pointLevel[f[fp]] <= faceLevel) {
          nAnchors++;
          if (nAnchors > 4) {
            protectedFace[facei] = true;
            break;
          }
        }
      }
    }

    syncTools::syncFaceList(*this, protectedFace, orEqOp<bool>());

    for (label facei = 0; facei < nInternalFaces(); facei++) {
      if (protectedFace[facei]) {
        protectedCell_.set(faceOwner()[facei], 1);
        nProtected++;
        protectedCell_.set(faceNeighbour()[facei], 1);
        nProtected++;
      }
    }

    for (label facei = nInternalFaces(); facei < nFaces(); facei++) {
      if (protectedFace[facei]) {
        protectedCell_.set(faceOwner()[facei], 1);
        nProtected++;
      }
    }

    // Also protect any cells that are less than hex
    forAll(cells(), celli) {
      const cell& cFaces = cells()[celli];

      if (cFaces.size() < 6) {
        if (protectedCell_.set(celli, 1)) nProtected++;
      } else {
        forAll(cFaces, cFacei) {
          if (faces()[cFaces[cFacei]].size() < 4) {
            if (protectedCell_.set(celli, 1)) nProtected++;
            break;
          }
        }
      }
    }

    // Check cells for 8 corner points
    checkEightAnchorPoints(protectedCell_, nProtected);
  }

  if (returnReduce(nProtected, sumOp<label>()) == 0)
    protectedCell_.clear();
  else {
    cellSet protectedCells(*this, "protectedCells", nProtected);
    forAll(protectedCell_, celli) if (protectedCell_[celli])
        protectedCells.insert(celli);

    Info << "Detected " << returnReduce(nProtected, sumOp<label>())
         << " cells that are protected from refinement."
         << " Writing these to cellSet " << protectedCells.name() << "."
         << endl;

    protectedCells.write();
  }
}

bool AMRFoam::updateAMR(const int& refineInterval, const int& maxRefinement,
                        volScalarField& vFld, const double& lowerRefineLevel,
                        const double& upperRefineLevel,
                        const double& unrefineAbove,
                        const double& unrefineBelow, const int& nBufferLayers,
                        const int &maxCells) {
  bool hasChanged = false;

  if (refineInterval == 0) {
    polyMesh::topoChanging(hasChanged);
    return hasChanged;
  } else if (refineInterval < 0) {
    FatalErrorInFunction << "Illegal refineInterval " << refineInterval << nl
                         << "The refineInterval value should"
                         << " be >= 1." << nl << exit(FatalError);
  }

  // Cannot refine at time = 0;
  if (fvMesh::time().timeIndex() > 0 &&
      fvMesh::time().timeIndex() % refineInterval == 0) {
    if (maxCells <= 0)
      FatalErrorInFunction << "Illegal maximum number of cells " << maxCells
                           << nl << "The maxCells value should"
                           << " be > 0." << nl << exit(FatalError);

    if (maxRefinement <= 0)
      FatalErrorInFunction << "Illegal maximum refinement level "
                           << maxRefinement << nl << "The maxCells value should"
                           << " be > 0." << nl << exit(FatalError);

    // Cells marked for refinement or otherwise protected from unrefinement.
    PackedBoolList refineCell(primitiveMesh::nCells());

    // Selecting refinement candidates based on lower and upper refinement
    // levels in scalar field
    selectRefineCandidates(lowerRefineLevel, upperRefineLevel, vFld,
                           refineCell);

    // Refinement Procedure Loop
    if (polyMesh::globalData().nTotalCells() < maxCells) {
      // Select subset of candidates. Take into account max allowable
      // cells, refinement level, protected cells.
      Foam::labelList cellsToRefine(
          selectRefineCells(maxCells, maxRefinement, refineCell));

      label nCellsToRefine = returnReduce(cellsToRefine.size(), sumOp<label>());

      if (nCellsToRefine > 0) {
        // Refine/update mesh and map fields
        autoPtr<mapPolyMesh> map = refineAMR(cellsToRefine);

        // Update refineCell. Note that some of the marked ones have not been
        // refined due to constraints.
        const labelList& cellMap = map().cellMap();
        const labelList& reverseCellMap = map().reverseCellMap();

        PackedBoolList newRefineCell(cellMap.size());

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
        for (label i = 0; i < nBufferLayers; i++) extendMarkedCells(refineCell);

        hasChanged = true;
      }
    }

    // Select unrefinable points that are not marked in refineCell
    labelList pointsToUnrefine(selectPts2Unrefine(
        unrefineAbove, unrefineBelow, refineCell, maxCellField(vFld)));

    label nSplitPoints = returnReduce(pointsToUnrefine.size(), sumOp<label>());

    if (nSplitPoints > 0) {
      // Refine/update mesh
      unrefineAMR(pointsToUnrefine);
      hasChanged = true;
    }

    if ((nRefinementIterations_ % 10) == 0) {
      // Unrefinement causes holes in the refinementHistory
      const_cast<refinementHistory&>(meshCutter().history()).compact();
    }
    nRefinementIterations_++;
  }

  polyMesh::topoChanging(hasChanged);
  if (hasChanged) polyMesh::moving(false);

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
    PackedBoolList refineCell(primitiveMesh::nCells());

    // Selecting refinement candidates based on lower and upper refinement
    // levels in scalar field
    selectRefineCandidates(0.5, 1.5, vFld, refineCell);

    // Refinement Procedure Loop
    if (polyMesh::globalData().nTotalCells() < maxCells) {
      // Select subset of candidates. Take into account max allowable
      // cells, refinement level, protected cells.
      Foam::labelList cellsToRefine(
          selectRefineCells(maxCells, maxRefinement, refineCell));

      label nCellsToRefine = returnReduce(cellsToRefine.size(), sumOp<label>());

      if (nCellsToRefine > 0) {
        // Refine/update mesh and map fields
        autoPtr<mapPolyMesh> map = refineAMR(cellsToRefine);

        // Update refineCell. Note that some of the marked ones have not been
        // refined due to constraints.
        const labelList& cellMap = map().cellMap();
        const labelList& reverseCellMap = map().reverseCellMap();

        PackedBoolList newRefineCell(cellMap.size());

        forAll(cellMap, celli) {
          label oldCelli = cellMap[celli];

          if (oldCelli < 0) {
            newRefineCell.set(celli, 1);
          } else if (reverseCellMap[oldCelli] != celli) {
            newRefineCell.set(celli, 1);
          } else {
            newRefineCell.set(celli, refineCell.get(oldCelli));
          }
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
        selectPts2Unrefine(-1, 0.5, refineCell, maxCellField(vFld)));

    label nSplitPoints = returnReduce(pointsToUnrefine.size(), sumOp<label>());
    autoPtr<mapPolyMesh> mapUnrefine;

    if (nSplitPoints > 0) {
      // Refine/update mesh
      mapUnrefine = unrefineAMR(pointsToUnrefine);

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

Foam::labelList AMRFoam::selectPts2Unrefine(
    const scalar unrefineAbove, const scalar unrefineBelow,
    const Foam::PackedBoolList& markedCell, const scalarField& pFld) {
  // All points that can be unrefined
  const labelList splitPoints(meshCutter_.getSplitPoints());

  DynamicList<label> newSplitPoints(splitPoints.size());

  forAll(splitPoints, i) {
    label pointi = splitPoints[i];

    if (pFld[pointi] < unrefineBelow || pFld[pointi] > unrefineAbove) {
      // Check that all cells are not marked
      const labelList& pCells = pointCells()[pointi];
      bool hasMarked = false;

      forAll(pCells, pCelli) {
        if (markedCell.get(pCells[pCelli])) {
          hasMarked = true;
          break;
        }
      }

      if (!hasMarked) newSplitPoints.append(pointi);
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

Foam::scalarField AMRFoam::cell2Pt(const Foam::scalarField& vFld) {
  Foam::scalarField pFld(nPoints());

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
  dynamicFvMesh::writeObject(time().writeFormat(), IOstream::currentVersion,
                             time().writeCompression(), true);
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

    for (int i = 0; i < myLines.size(); i++)
      incomingField[i] = std::atof(nemAux::strToChar(myLines[i]).get());

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
  for (int i = 0; i < vec.size(); i++) vFld[i] = vec[i];
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

    for (int i = 0; i < myLines.size(); i++)
      incomingField[i] = std::atof(nemAux::strToChar(myLines[i]).get());

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

  returnGrad = Foam::fvc::grad(fldGrd);
  returnGrad.write();
  volScalarField returnScalar = initGradientField("gradient");
  returnScalar = Foam::mag(returnGrad);

  return returnScalar;
}

Foam::label AMRFoam::count(const PackedBoolList& l, const unsigned int val) {
  label n = 0;
  forAll(l, i) if (l.get(i) == val) n++;
  return n;
}

void AMRFoam::calculateProtectedCells(PackedBoolList& unrefineableCell) const {
  if (protectedCell_.empty()) {
    unrefineableCell.clear();
    return;
  }

  const labelList& cellLevel = meshCutter_.cellLevel();
  unrefineableCell = protectedCell_;

  // Get neighbouring cell level
  labelList neiLevel(nFaces() - nInternalFaces());

  for (label facei = nInternalFaces(); facei < nFaces(); facei++)
    neiLevel[facei - nInternalFaces()] = cellLevel[faceOwner()[facei]];
  syncTools::swapBoundaryFaceList(*this, neiLevel);

  while (true) {
    // Pick up faces on border of protected cells
    boolList seedFace(nFaces(), false);
    forAll(faceNeighbour(), facei) {
      label own = faceOwner()[facei];
      bool ownProtected = unrefineableCell.get(own);
      label nei = faceNeighbour()[facei];
      bool neiProtected = unrefineableCell.get(nei);

      if (ownProtected && (cellLevel[nei] > cellLevel[own]))
        seedFace[facei] = true;
      else if (neiProtected && (cellLevel[own] > cellLevel[nei]))
        seedFace[facei] = true;
    }
    for (label facei = nInternalFaces(); facei < nFaces(); facei++) {
      label own = faceOwner()[facei];
      bool ownProtected = unrefineableCell.get(own);
      if (ownProtected && (neiLevel[facei - nInternalFaces()] > cellLevel[own]))
        seedFace[facei] = true;
    }

    syncTools::syncFaceList(*this, seedFace, orEqOp<bool>());

    // Extend unrefineableCell
    bool hasExtended = false;

    for (label facei = 0; facei < nInternalFaces(); facei++) {
      if (seedFace[facei]) {
        label own = faceOwner()[facei];
        if (unrefineableCell.get(own) == 0) {
          unrefineableCell.set(own, 1);
          hasExtended = true;
        }
        label nei = faceNeighbour()[facei];
        if (unrefineableCell.get(nei) == 0) {
          unrefineableCell.set(nei, 1);
          hasExtended = true;
        }
      }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); facei++) {
      if (seedFace[facei]) {
        label own = faceOwner()[facei];
        if (unrefineableCell.get(own) == 0) {
          unrefineableCell.set(own, 1);
          hasExtended = true;
        }
      }
    }

    if (!returnReduce(hasExtended, orOp<bool>())) break;
  }
}

void AMRFoam::readDict() {
  dictionary refineDict(
      IOdictionary(IOobject("dynamicMeshDict", time().constant(), *this,
                            IOobject::MUST_READ_IF_MODIFIED, IOobject::NO_WRITE,
                            false))
          .optionalSubDict(typeName + "Coeffs"));

  List<Pair<word>> fluxVelocities =
      List<Pair<word>>(refineDict.lookup("correctFluxes"));
  // Rework into hashtable.
  correctFluxes_.resize(fluxVelocities.size());
  forAll(fluxVelocities, i)
      correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);

  dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
}

// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh> Foam::AMRFoam::refineAMR(
    const labelList& cellsToRefine) {
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
    for (label facei = 0; facei < nInternalFaces(); facei++) {
      label oldFacei = map().faceMap()[facei];
      if (oldFacei >= nInternalFaces())
        FatalErrorInFunction << "New internal face:" << facei
                             << " fc:" << faceCentres()[facei]
                             << " originates from boundary oldFace:" << oldFacei
                             << abort(FatalError);
    }
  }

  // Update fields
  updateMesh(map);

  // Move mesh
  // pointField newPoints;
  // if (map().hasMotionPoints())
  // {
  //     newPoints = map().preMotionPoints();
  // }
  // else
  // {
  //     newPoints = points();
  // }
  // movePoints(newPoints);
  // motionPtr_->updateMesh(map);
  // movePoints(motionPtr_->newPoints());

  // Correct the flux for modified/added faces. All the faces which only
  // have been renumbered will already have been handled by the mapping.
  {
    const labelList& faceMap = map().faceMap();
    const labelList& reverseFaceMap = map().reverseFaceMap();

    // Storage for any master faces. These will be the original faces
    // on the coarse cell that get split into four (or rather the
    // master face gets modified and three faces get added from the master)
    labelHashSet masterFaces(4 * cellsToRefine.size());

    forAll(faceMap, facei) {
      label oldFacei = faceMap[facei];

      if (oldFacei >= 0) {
        label masterFacei = reverseFaceMap[oldFacei];

        if (masterFacei < 0)
          FatalErrorInFunction << "Problem: should not have removed faces"
                               << " when refining." << nl << "face:" << facei
                               << abort(FatalError);
        else if (masterFacei != facei)
          masterFaces.insert(masterFacei);
      }
    }

    if (debug)
      Pout << "Found " << masterFaces.size() << " split faces " << endl;

    HashTable<surfaceScalarField*> fluxes(lookupClass<surfaceScalarField>());
    forAllIter(HashTable<surfaceScalarField*>, fluxes, iter) {
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

      // const word& UName = correctFluxes_[iter.key()];
      const word& UName = "U";
      if (UName == "none") continue;

      if (UName == "NaN") {
        Pout << "Setting surfaceScalarField " << iter.key() << " to NaN"
             << endl;

        surfaceScalarField& phi = *iter();
        sigFpe::fillNan(phi.primitiveFieldRef());

        continue;
      }

      if (debug)
        Pout << "Mapping flux " << iter.key() << " using interpolated flux "
             << UName << endl;

      surfaceScalarField& phi = *iter();
      const surfaceScalarField phiU(
          fvc::interpolate(lookupObject<volVectorField>(UName)) & Sf());

      // Recalculate new internal faces.
      for (label facei = 0; facei < nInternalFaces(); facei++) {
        label oldFacei = faceMap[facei];

        if (oldFacei == -1)
          phi[facei] = phiU[facei];  // Inflated/appended
        else if (reverseFaceMap[oldFacei] != facei)
          phi[facei] = phiU[facei];  // face-from-masterface
      }

      // Recalculate new boundary faces.
      surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

      forAll(phiBf, patchi) {
        fvsPatchScalarField& patchPhi = phiBf[patchi];
        const fvsPatchScalarField& patchPhiU = phiU.boundaryField()[patchi];
        label facei = patchPhi.patch().start();

        forAll(patchPhi, i) {
          label oldFacei = faceMap[facei];

          if (oldFacei == -1)
            patchPhi[i] = patchPhiU[i];  // Inflated/appended
          else if (reverseFaceMap[oldFacei] != facei)
            patchPhi[i] = patchPhiU[i];  // face-from-masterface

          facei++;
        }
      }

      // Update master faces
      forAllConstIter(labelHashSet, masterFaces, iter) {
        label facei = iter.key();

        if (isInternalFace(facei))
          phi[facei] = phiU[facei];
        else {
          label patchi = boundaryMesh().whichPatch(facei);
          label i = facei - boundaryMesh()[patchi].start();
          const fvsPatchScalarField& patchPhiU = phiU.boundaryField()[patchi];
          fvsPatchScalarField& patchPhi = phiBf[patchi];
          patchPhi[i] = patchPhiU[i];
        }
      }
    }
  }

  // Update numbering of cells/vertices.
  meshCutter_.updateMesh(map);

  // Update numbering of protectedCell_
  if (protectedCell_.size()) {
    PackedBoolList newProtectedCell(nCells());
    forAll(newProtectedCell, celli) {
      label oldCelli = map().cellMap()[celli];
      newProtectedCell.set(celli, protectedCell_.get(oldCelli));
    }
    protectedCell_.transfer(newProtectedCell);
  }

  // Debug: Check refinement levels (across faces only)
  meshCutter_.checkRefinementLevels(-1, labelList(0));

  return map;
}

Foam::autoPtr<Foam::mapPolyMesh> Foam::AMRFoam::unrefineAMR(
    const labelList& splitPoints) {
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
    forAll(splitPoints, i) {
      label pointi = splitPoints[i];
      const labelList& pEdges = pointEdges()[pointi];

      forAll(pEdges, j) {
        label otherPointi = edges()[pEdges[j]].otherVertex(pointi);
        const labelList& pFaces = pointFaces()[otherPointi];

        forAll(pFaces, pFacei)
            faceToSplitPoint.insert(pFaces[pFacei], otherPointi);
      }
    }
  }

  // Change mesh and generate map.
  // autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
  autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

  Info << "Unrefined from " << returnReduce(map().nOldCells(), sumOp<label>())
       << " to " << globalData().nTotalCells() << " cells." << endl;

  // Update fields
  updateMesh(map);

  // Move mesh
  // pointField newPoints;
  // if (map().hasMotionPoints())
  // {
  //     newPoints = map().preMotionPoints();
  // }
  // else
  // {
  //     newPoints = points();
  // }
  // movePoints(newPoints);
  // motionPtr_->updateMesh(map);
  // movePoints(motionPtr_->newPoints());

  // Correct the flux for modified faces.
  {
    const labelList& reversePointMap = map().reversePointMap();
    const labelList& reverseFaceMap = map().reverseFaceMap();

    HashTable<surfaceScalarField*> fluxes(lookupClass<surfaceScalarField>());
    forAllIter(HashTable<surfaceScalarField*>, fluxes, iter) {
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

      // const word& UName = correctFluxes_[iter.key()];
      const word& UName = "U";

      if (UName == "none") continue;

      if (debug)
        Info << "Mapping flux " << iter.key() << " using interpolated flux "
             << UName << endl;

      surfaceScalarField& phi = *iter();
      surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

      const surfaceScalarField phiU(
          fvc::interpolate(lookupObject<volVectorField>(UName)) & Sf());

      forAllConstIter(Map<label>, faceToSplitPoint, iter) {
        label oldFacei = iter.key();
        label oldPointi = iter();

        if (reversePointMap[oldPointi] < 0) {
          // midpoint was removed. See if face still exists.
          label facei = reverseFaceMap[oldFacei];

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
  meshCutter_.updateMesh(map);

  // Update numbering of protectedCell_
  if (protectedCell_.size()) {
    PackedBoolList newProtectedCell(nCells());

    forAll(newProtectedCell, celli) {
      label oldCelli = map().cellMap()[celli];
      if (oldCelli >= 0)
        newProtectedCell.set(celli, protectedCell_.get(oldCelli));
    }

    protectedCell_.transfer(newProtectedCell);
  }

  // Debug: Check refinement levels (across faces only)
  meshCutter_.checkRefinementLevels(-1, labelList(0));

  return map;
}

Foam::scalarField AMRFoam::maxPointField(const scalarField& pFld) const {
  scalarField vFld(nCells(), -GREAT);
  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];
    forAll(pCells, i) vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointi]);
  }
  return vFld;
}

Foam::scalarField AMRFoam::maxCellField(const volScalarField& vFld) const {
  scalarField pFld(nPoints(), -GREAT);
  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];
    forAll(pCells, i) pFld[pointi] = max(pFld[pointi], vFld[pCells[i]]);
  }
  return pFld;
}

Foam::scalarField AMRFoam::cellToPoint(const scalarField& vFld) const {
  scalarField pFld(nPoints());

  forAll(pointCells(), pointi) {
    const labelList& pCells = pointCells()[pointi];
    scalar sum = 0.0;
    forAll(pCells, i) sum += vFld[pCells[i]];

    pFld[pointi] = sum / pCells.size();
  }
  return pFld;
}

Foam::scalarField AMRFoam::error(const scalarField& fld, const scalar minLevel,
                                 const scalar maxLevel) const {
  scalarField c(fld.size(), -1);
  forAll(fld, i) {
    scalar err = min(fld[i] - minLevel, maxLevel - fld[i]);
    if (err >= 0) c[i] = err;
  }
  return c;
}

void AMRFoam::selectRefineCandidates(const scalar lowerRefineLevel,
                                     const scalar upperRefineLevel,
                                     const scalarField& vFld,
                                     PackedBoolList& candidateCell) const {
  // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
  // higher more desirable to be refined).
  scalarField cellError(maxPointField(
      error(cellToPoint(vFld), lowerRefineLevel, upperRefineLevel)));

  // Mark cells that are candidates for refinement.
  forAll(cellError, celli) if (cellError[celli] > 0)
      candidateCell.set(celli, 1);
}

Foam::labelList AMRFoam::selectRefineCells(
    const label maxCells, const label maxRefinement,
    const PackedBoolList& candidateCell) const {
  // Every refined cell causes 7 extra cells
  label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

  const labelList& cellLevel = meshCutter_.cellLevel();

  // Mark cells that cannot be refined since they would trigger refinement
  // of protected cells (since 2:1 cascade)
  PackedBoolList unrefineableCell;
  calculateProtectedCells(unrefineableCell);

  // Count current selection
  label nLocalCandidates = count(candidateCell, 1);
  label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

  // Collect all cells
  DynamicList<label> candidates(nLocalCandidates);

  if (nCandidates < nTotToRefine) {
    forAll(candidateCell, celli) {
      if (cellLevel[celli] < maxRefinement && candidateCell.get(celli) &&
          (unrefineableCell.empty() || !unrefineableCell.get(celli))) {
        candidates.append(celli);
      }
    }
  } else {
    // Sort by error? For now just truncate.
    for (label level = 0; level < maxRefinement; level++) {
      forAll(candidateCell, celli) {
        if (cellLevel[celli] == level && candidateCell.get(celli) &&
            (unrefineableCell.empty() || !unrefineableCell.get(celli))) {
          candidates.append(celli);
        }
      }

      if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine) break;
    }
  }

  // Guarantee 2:1 refinement after refinement
  labelList consistentSet(
      meshCutter_.consistentRefinement(candidates.shrink(), true));

  Info << "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
       << " cells for refinement out of " << globalData().nTotalCells() << "."
       << endl;

  return consistentSet;
}

void AMRFoam::extendMarkedCells(PackedBoolList& markedCell) const {
  // Mark faces using any marked cell
  boolList markedFace(nFaces(), false);

  forAll(markedCell, celli) {
    if (markedCell.get(celli)) {
      const cell& cFaces = cells()[celli];
      forAll(cFaces, i) markedFace[cFaces[i]] = true;
    }
  }

  syncTools::syncFaceList(*this, markedFace, orEqOp<bool>());

  // Update cells using any markedFace
  for (label facei = 0; facei < nInternalFaces(); facei++) {
    if (markedFace[facei]) {
      markedCell.set(faceOwner()[facei], 1);
      markedCell.set(faceNeighbour()[facei], 1);
    }
  }
  for (label facei = nInternalFaces(); facei < nFaces(); facei++) {
    if (markedFace[facei]) {
      markedCell.set(faceOwner()[facei], 1);
    }
  }
}

void AMRFoam::checkEightAnchorPoints(PackedBoolList& protectedCell,
                                     label& nProtected) const {
  const labelList& cellLevel = meshCutter_.cellLevel();
  const labelList& pointLevel = meshCutter_.pointLevel();

  labelList nAnchorPoints(nCells(), 0);

  forAll(pointLevel, pointi) {
    const labelList& pCells = pointCells(pointi);

    forAll(pCells, pCelli) {
      label celli = pCells[pCelli];

      if (pointLevel[pointi] <= cellLevel[celli]) {
        // Check if cell has already 8 anchor points -> protect cell
        if (nAnchorPoints[celli] == 8)
          if (protectedCell.set(celli, true)) nProtected++;

        if (!protectedCell[celli]) nAnchorPoints[celli]++;
      }
    }
  }

  forAll(protectedCell, celli) {
    if (!protectedCell[celli] && nAnchorPoints[celli] != 8) {
      protectedCell.set(celli, true);
      nProtected++;
    }
  }
}

bool AMRFoam::writeObject(IOstream::streamFormat fmt,
                          IOstream::versionNumber ver,
                          IOstream::compressionType cmp,
                          const bool valid) const {
  // Force refinement data to go to the current time directory.
  const_cast<hexRef8&>(meshCutter_).setInstance(time().timeName());

  bool writeOk = (dynamicFvMesh::writeObject(fmt, ver, cmp, valid) &&
                  meshCutter_.write(valid));

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

}  // namespace Foam

// ************************************************************************* //
