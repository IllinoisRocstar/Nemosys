#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "createPatch.H"
#include "getDicts.H"

using namespace Foam;

createPatch::createPatch() {}

createPatch::~createPatch() {}

void createPatch::execute(const int& dirStat, const std::string& firstPtch,
                          const std::unique_ptr<Foam::dictionary>& cpDict) {
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);
  auto fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  auto fvSolution_ = initFoam->createFvSolution(writeDicts);

  // Create time class without reading controlDict
  Foam::Time runTime(controlDict_.get(), ".", ".");
  Foam::argList::noParallel();

  std::string secondPtch;

  if (dirStat == 1)
    secondPtch = "domain2";
  else
    secondPtch = "domain1";

  std::string regnName;

  // This loop wraps around utility execution to merge patches of all packs and
  // surounding regions into two distinct patches.
  for (int i = 0; i < 2; i++) {
    if (i == 0) {
      regnName = firstPtch;
    }
    if (i == 1) {
      regnName = secondPtch;
    }
    const bool overwrite = true;

    Foam::word meshRegionName = regnName;

    Foam::word regionName;
    regionName = regnName;
    Foam::Info << "Create polyMesh for time = " << runTime.timeName()
               << Foam::nl << Foam::endl;

    Foam::polyMesh mesh(Foam::IOobject(regionName, runTime.timeName(), runTime,
                                       Foam::IOobject::MUST_READ));

    const word oldInstance = mesh.pointsInstance();

    const word dictName("createPatchDict");

    IOobject dictIO(dictName, runTime.system(), mesh, IOobject::MUST_READ,
                    IOobject::NO_WRITE);

    Info << "Reading " << dictName << nl << endl;

    IOdictionary dict(dictIO);

    // Whether to synchronise points
    const Switch pointSync(dict.lookup("pointSync"));

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // If running parallel check same patches everywhere
    patches.checkParallelSync(true);

    dumpCyclicMatch("initial_", mesh);

    // Read patch construct info from dictionary
    PtrList<dictionary> patchSources(dict.lookup("patches"));

    HashSet<word> addedPatchNames;
    forAll(patchSources, addedI) {
      const dictionary& dict = patchSources[addedI];
      addedPatchNames.insert(dict.get<word>("name"));
    }

    // 1. Add all new patches
    // ~~~~~~~~~~~~~~~~~~~~~~

    if (patchSources.size()) {
      // Old and new patches.
      DynamicList<polyPatch*> allPatches(patches.size() + patchSources.size());

      label startFacei = mesh.nInternalFaces();

      // Copy old patches.
      forAll(patches, patchi) {
        const polyPatch& pp = patches[patchi];

        if (!isA<processorPolyPatch>(pp)) {
          allPatches.append(
              pp.clone(patches, patchi, pp.size(), startFacei).ptr());
          startFacei += pp.size();
        }
      }

      forAll(patchSources, addedI) {
        const dictionary& dict = patchSources[addedI];

        word patchName(dict.lookup("name"));

        label destPatchi = patches.findPatchID(patchName);

        if (destPatchi == -1) {
          dictionary patchDict(dict.subDict("patchInfo"));

          destPatchi = allPatches.size();

          Info << "Adding new patch " << patchName << " as patch " << destPatchi
               << " from " << patchDict << endl;

          patchDict.set("nFaces", 0);
          patchDict.set("startFace", startFacei);

          // Add an empty patch.
          allPatches.append(
              polyPatch::New(patchName, patchDict, destPatchi, patches).ptr());
        } else {
          Info << "Patch '" << patchName << "' already exists.  Only "
               << "moving patch faces - type will remain the same" << endl;
        }
      }

      // Copy old patches.
      forAll(patches, patchi) {
        const polyPatch& pp = patches[patchi];

        if (isA<processorPolyPatch>(pp)) {
          allPatches.append(
              pp.clone(patches, patchi, pp.size(), startFacei).ptr());
          startFacei += pp.size();
        }
      }

      allPatches.shrink();
      mesh.removeBoundary();
      mesh.addPatches(allPatches);
      Info << endl;
    }

    // 2. Repatch faces
    // ~~~~~~~~~~~~~~~~

    polyTopoChange meshMod(mesh);

    forAll(patchSources, addedI) {
      const dictionary& dict = patchSources[addedI];

      const word patchName(dict.lookup("name"));
      label destPatchi = patches.findPatchID(patchName);

      if (destPatchi == -1) {
        FatalErrorInFunction << "patch " << patchName << " not added. Problem."
                             << abort(FatalError);
      }

      const word sourceType(dict.lookup("constructFrom"));

      if (sourceType == "patches") {
        labelHashSet patchSources(
            patches.patchSet(wordReList(dict.lookup("patches"))));

        // Repatch faces of the patches.
        forAllConstIter(labelHashSet, patchSources, iter) {
          const polyPatch& pp = patches[iter.key()];

          Info << "Moving faces from patch " << pp.name() << " to patch "
               << destPatchi << endl;

          forAll(pp, i) {
            changePatchID(mesh, pp.start() + i, destPatchi, meshMod);
          }
        }
      } else if (sourceType == "set") {
        const word setName(dict.lookup("set"));

        faceSet faces(mesh, setName);

        Info << "Read " << returnReduce(faces.size(), sumOp<label>())
             << " faces from faceSet " << faces.name() << endl;

        // Sort (since faceSet contains faces in arbitrary order)
        labelList faceLabels(faces.toc());

        SortableList<label> patchFaces(faceLabels);

        forAll(patchFaces, i) {
          label facei = patchFaces[i];

          if (mesh.isInternalFace(facei)) {
            FatalErrorInFunction
                << "Face " << facei << " specified in set " << faces.name()
                << " is not an external face of the mesh." << endl
                << "This application can only repatch existing boundary"
                << " faces." << exit(FatalError);
          }

          changePatchID(mesh, facei, destPatchi, meshMod);
        }
      } else {
        FatalErrorInFunction << "Invalid source type " << sourceType << endl
                             << "Valid source types are 'patches' 'set'"
                             << exit(FatalError);
      }
    }
    Info << endl;

    // Change mesh, use inflation to reforce calculation of transformation
    // tensors.
    Info << "Doing topology modification to order faces." << nl << endl;
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, true);
    mesh.movePoints(map().preMotionPoints());

    dumpCyclicMatch("coupled_", mesh);

    // Synchronise points.
    if (!pointSync) {
      Info << "Not synchronising points." << nl << endl;
    } else {
      Info << "Synchronising points." << nl << endl;

      // This is a bit tricky. Both normal and position might be out and
      // current separation also includes the normal
      // ( separation_ = (nf&(Cr - Cf))*nf ).

      // For cyclic patches:
      // - for separated ones use user specified offset vector

      forAll(mesh.boundaryMesh(), patchi) {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if (pp.size() && isA<coupledPolyPatch>(pp)) {
          const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>(pp);

          if (cpp.separated()) {
            Info << "On coupled patch " << pp.name() << " separation[0] was "
                 << cpp.separation()[0] << endl;

            if (isA<cyclicPolyPatch>(pp) && pp.size()) {
              const cyclicPolyPatch& cycpp = refCast<const cyclicPolyPatch>(pp);

              if (cycpp.transform() == cyclicPolyPatch::TRANSLATIONAL) {
                // Force to wanted separation
                Info << "On cyclic translation patch " << pp.name()
                     << " forcing uniform separation of "
                     << cycpp.separationVector() << endl;
                const_cast<vectorField&>(cpp.separation()) =
                    pointField(1, cycpp.separationVector());
              } else {
                const cyclicPolyPatch& nbr = cycpp.neighbPatch();
                const_cast<vectorField&>(cpp.separation()) =
                    pointField(1, nbr[0].centre(mesh.points()) -
                                      cycpp[0].centre(mesh.points()));
              }
            }
            Info << "On coupled patch " << pp.name()
                 << " forcing uniform separation of " << cpp.separation()
                 << endl;
          } else if (!cpp.parallel()) {
            Info << "On coupled patch " << pp.name()
                 << " forcing uniform rotation of " << cpp.forwardT()[0]
                 << endl;

            const_cast<tensorField&>(cpp.forwardT()).setSize(1);
            const_cast<tensorField&>(cpp.reverseT()).setSize(1);

            Info << "On coupled patch " << pp.name()
                 << " forcing uniform rotation of " << cpp.forwardT() << endl;
          }
        }
      }

      Info << "Synchronising points." << endl;

      pointField newPoints(mesh.points());

      syncPoints(mesh, newPoints, minMagSqrEqOp<vector>(),
                 point(GREAT, GREAT, GREAT));

      scalarField diff(mag(newPoints - mesh.points()));
      Info << "Points changed by average:" << gAverage(diff)
           << " max:" << gMax(diff) << nl << endl;

      mesh.movePoints(newPoints);
    }

    // 3. Remove zeros-sized patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Info << "Removing patches with no faces in them." << nl << endl;
    filterPatches(mesh, addedPatchNames);

    dumpCyclicMatch("final_", mesh);

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    if (!overwrite) {
      runTime++;
    } else {
      mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info << "Writing repatched mesh to " << runTime.timeName() << nl << endl;
    mesh.write();
  }

  Info << "End\n" << endl;
}

// createPatch
void createPatch::changePatchID(const Foam::polyMesh& mesh,
                                const Foam::label faceID,
                                const Foam::label patchID,
                                Foam::polyTopoChange& meshMod) {
  const label zoneID = mesh.faceZones().whichZone(faceID);

  bool zoneFlip = false;

  if (zoneID >= 0) {
    const faceZone& fZone = mesh.faceZones()[zoneID];

    zoneFlip = fZone.flipMap()[fZone.whichFace(faceID)];
  }

  meshMod.setAction(polyModifyFace(mesh.faces()[faceID],      // face
                                   faceID,                    // face ID
                                   mesh.faceOwner()[faceID],  // owner
                                   -1,                        // neighbour
                                   false,                     // flip flux
                                   patchID,                   // patch ID
                                   false,    // remove from zone
                                   zoneID,   // zone ID
                                   zoneFlip  // zone flip
                                   ));
}

void createPatch::filterPatches(
    Foam::polyMesh& mesh, const Foam::HashSet<Foam::word>& addedPatchNames) {
  const polyBoundaryMesh& patches = mesh.boundaryMesh();

  // Patches to keep
  DynamicList<polyPatch*> allPatches(patches.size());

  label nOldPatches = returnReduce(patches.size(), sumOp<label>());

  // Copy old patches.
  forAll(patches, patchi) {
    const polyPatch& pp = patches[patchi];

    // Note: reduce possible since non-proc patches guaranteed in same order
    if (!isA<processorPolyPatch>(pp)) {
      // Add if
      // - non zero size
      // - or added from the createPatchDict
      // - or cyclic (since referred to by other cyclic half or
      //   proccyclic)

      if (addedPatchNames.found(pp.name()) ||
          returnReduce(pp.size(), sumOp<label>()) > 0 ||
          isA<coupledPolyPatch>(pp)) {
        allPatches.append(
            pp.clone(patches, allPatches.size(), pp.size(), pp.start()).ptr());
      } else {
        Info << "Removing zero-sized patch " << pp.name() << " type "
             << pp.type() << " at position " << patchi << endl;
      }
    }
  }
  // Copy non-empty processor patches
  forAll(patches, patchi) {
    const polyPatch& pp = patches[patchi];

    if (isA<processorPolyPatch>(pp)) {
      if (pp.size()) {
        allPatches.append(
            pp.clone(patches, allPatches.size(), pp.size(), pp.start()).ptr());
      } else {
        Info << "Removing empty processor patch " << pp.name()
             << " at position " << patchi << endl;
      }
    }
  }

  label nAllPatches = returnReduce(allPatches.size(), sumOp<label>());
  if (nAllPatches != nOldPatches) {
    Info << "Removing patches." << endl;
    allPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(allPatches);
  } else {
    Info << "No patches removed." << endl;
    forAll(allPatches, i) { delete allPatches[i]; }
  }
}

void createPatch::dumpCyclicMatch(const Foam::fileName& prefix,
                                  const Foam::polyMesh& mesh) {
  const polyBoundaryMesh& patches = mesh.boundaryMesh();

  forAll(patches, patchi) {
    if (isA<cyclicPolyPatch>(patches[patchi]) &&
        refCast<const cyclicPolyPatch>(patches[patchi]).owner()) {
      const cyclicPolyPatch& cycPatch =
          refCast<const cyclicPolyPatch>(patches[patchi]);

      // Dump patches
      {
        OFstream str(prefix + cycPatch.name() + ".obj");
        Pout << "Dumping " << cycPatch.name() << " faces to " << str.name()
             << endl;
        meshTools::writeOBJ(str, cycPatch, cycPatch.points());
      }

      const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
      {
        OFstream str(prefix + nbrPatch.name() + ".obj");
        Pout << "Dumping " << nbrPatch.name() << " faces to " << str.name()
             << endl;
        meshTools::writeOBJ(str, nbrPatch, nbrPatch.points());
      }

      // Lines between corresponding face centres
      OFstream str(prefix + cycPatch.name() + nbrPatch.name() + "_match.obj");
      label vertI = 0;

      Pout << "Dumping cyclic match as lines between face centres to "
           << str.name() << endl;

      forAll(cycPatch, facei) {
        const point& fc0 = mesh.faceCentres()[cycPatch.start() + facei];
        meshTools::writeOBJ(str, fc0);
        vertI++;
        const point& fc1 = mesh.faceCentres()[nbrPatch.start() + facei];
        meshTools::writeOBJ(str, fc1);
        vertI++;

        str << "l " << vertI - 1 << ' ' << vertI << nl;
      }
    }
  }
}

void createPatch::separateList(const Foam::vectorField& separation,
                               Foam::UList<Foam::vector>& field) {
  if (separation.size() == 1) {
    // Single value for all.

    forAll(field, i) { field[i] += separation[0]; }
  } else if (separation.size() == field.size()) {
    forAll(field, i) { field[i] += separation[i]; }
  } else {
    FatalErrorInFunction
        << "Sizes of field and transformation not equal. field:" << field.size()
        << " transformation:" << separation.size() << abort(FatalError);
  }
}

template <class CombineOp>
void createPatch::syncPoints(const Foam::polyMesh& mesh,
                             Foam::pointField& points, const CombineOp& cop,
                             const Foam::point& nullValue) {
  if (points.size() != mesh.nPoints()) {
    FatalErrorInFunction << "Number of values " << points.size()
                         << " is not equal to the number of points in the mesh "
                         << mesh.nPoints() << abort(FatalError);
  }

  const polyBoundaryMesh& patches = mesh.boundaryMesh();

  // Is there any coupled patch with transformation?
  bool hasTransformation = false;

  if (Pstream::parRun()) {
    // Send

    forAll(patches, patchi) {
      const polyPatch& pp = patches[patchi];

      if (isA<processorPolyPatch>(pp) && pp.nPoints() > 0 &&
          refCast<const processorPolyPatch>(pp).owner()) {
        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(pp);

        // Get data per patchPoint in neighbouring point numbers.
        pointField patchInfo(procPatch.nPoints(), nullValue);

        const labelList& meshPts = procPatch.meshPoints();
        const labelList& nbrPts = procPatch.neighbPoints();

        forAll(nbrPts, pointi) {
          label nbrPointi = nbrPts[pointi];
          if (nbrPointi >= 0 && nbrPointi < patchInfo.size()) {
            patchInfo[nbrPointi] = points[meshPts[pointi]];
          }
        }

        OPstream toNbr(Pstream::commsTypes::blocking, procPatch.neighbProcNo());
        toNbr << patchInfo;
      }
    }

    // Receive and set.

    forAll(patches, patchi) {
      const polyPatch& pp = patches[patchi];

      if (isA<processorPolyPatch>(pp) && pp.nPoints() > 0 &&
          !refCast<const processorPolyPatch>(pp).owner()) {
        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(pp);

        pointField nbrPatchInfo(procPatch.nPoints());
        {
          // We do not know the number of points on the other side
          // so cannot use Pstream::read.
          IPstream fromNbr(Pstream::commsTypes::blocking,
                           procPatch.neighbProcNo());
          fromNbr >> nbrPatchInfo;
        }
        // Null any value which is not on neighbouring processor
        nbrPatchInfo.setSize(procPatch.nPoints(), nullValue);

        if (!procPatch.parallel()) {
          hasTransformation = true;
          transformList(procPatch.forwardT(), nbrPatchInfo);
        } else if (procPatch.separated()) {
          hasTransformation = true;
          separateList(-procPatch.separation(), nbrPatchInfo);
        }

        const labelList& meshPts = procPatch.meshPoints();

        forAll(meshPts, pointi) {
          label meshPointi = meshPts[pointi];
          points[meshPointi] = nbrPatchInfo[pointi];
        }
      }
    }
  }

  // Do the cyclics.
  forAll(patches, patchi) {
    const polyPatch& pp = patches[patchi];

    if (isA<cyclicPolyPatch>(pp) &&
        refCast<const cyclicPolyPatch>(pp).owner()) {
      const cyclicPolyPatch& cycPatch = refCast<const cyclicPolyPatch>(pp);

      const edgeList& coupledPoints = cycPatch.coupledPoints();
      const labelList& meshPts = cycPatch.meshPoints();
      const cyclicPolyPatch& nbrPatch = cycPatch.neighbPatch();
      const labelList& nbrMeshPts = nbrPatch.meshPoints();

      pointField half0Values(coupledPoints.size());

      forAll(coupledPoints, i) {
        const edge& e = coupledPoints[i];
        label point0 = meshPts[e[0]];
        half0Values[i] = points[point0];
      }

      if (!cycPatch.parallel()) {
        hasTransformation = true;
        transformList(cycPatch.reverseT(), half0Values);
      } else if (cycPatch.separated()) {
        hasTransformation = true;
        separateList(cycPatch.separation(), half0Values);
      }

      forAll(coupledPoints, i) {
        const edge& e = coupledPoints[i];
        label point1 = nbrMeshPts[e[1]];
        points[point1] = half0Values[i];
      }
    }
  }

  //- Note: hasTransformation is only used for warning messages so
  //  reduction not strictly nessecary.
  // reduce(hasTransformation, orOp<bool>());

  // Synchronize multiple shared points.
  const globalMeshData& pd = mesh.globalData();

  if (pd.nGlobalPoints() > 0) {
    if (hasTransformation) {
      WarningInFunction << "There are decomposed cyclics in this mesh with"
                        << " transformations." << endl
                        << "This is not supported. The result will be incorrect"
                        << endl;
    }

    // Values on shared points.
    pointField sharedPts(pd.nGlobalPoints(), nullValue);

    forAll(pd.sharedPointLabels(), i) {
      label meshPointi = pd.sharedPointLabels()[i];
      // Fill my entries in the shared points
      sharedPts[pd.sharedPointAddr()[i]] = points[meshPointi];
    }

    // Combine on master.
    Pstream::listCombineGather(sharedPts, cop);
    Pstream::listCombineScatter(sharedPts);

    // Now we will all have the same information. Merge it back with
    // my local information.
    forAll(pd.sharedPointLabels(), i) {
      label meshPointi = pd.sharedPointLabels()[i];
      points[meshPointi] = sharedPts[pd.sharedPointAddr()[i]];
    }
  }
}
