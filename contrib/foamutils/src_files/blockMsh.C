#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include "blockMsh.H"
#include "getDicts.H"

using namespace Foam;
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

blockMsh::blockMsh() {}

blockMsh::~blockMsh() {}

// Implementation of blockMesh code
Foam::fvMesh* blockMsh::generate(
    const std::unique_ptr<Foam::Time>& runTime,
    const std::unique_ptr<Foam::dictionary>& blockMshDict_, const bool& write) {
  if (!runTime) {
    std::cerr << "Cannot find OpenFOAM runtime, exiting ...." << std::endl;
    throw;
  }
  if (!blockMshDict_) {
    std::cerr << "Cannot find blockMeshDict, exiting ...." << std::endl;
    throw;
  }
  Foam::argList::noParallel();

  word regionName;
  word regionPath;
  regionName = polyMesh::defaultRegion;
  regionPath = regionName;

  // Locating blockMeshDict in system directory
  const word dictName("blockMeshDict");
  fileName dictPath;
  dictPath = runTime->system() / dictName;

  // Looks for and cleans polyMesh (Can use boost directory if needed)
  fileName polyMeshPath(runTime->path() / runTime->constant() /
                        polyMesh::meshSubDir);

  if (exists(polyMeshPath)) {
    if (exists(polyMeshPath / dictName)) {
      Info << "Not deleting polyMesh directory " << nl << "    " << polyMeshPath
           << nl << "    because it contains " << dictName << endl;
    } else {
      Info << "Deleting polyMesh directory" << nl << "    " << polyMeshPath
           << endl;
      rmDir(polyMeshPath);
    }
  }

  // Creates blockMesh from defined dictionary. Change it to NO_READ for IN
  // MEMORY
  IOobject meshDictIO(dictPath, *runTime, IOobject::NO_READ, IOobject::NO_WRITE,
                      true);

  IOdictionary meshDict(meshDictIO, blockMshDict_.get());  // For IN MEMORY
  // IOdictionary meshDict(meshDictIO);

  Info << "Creating block mesh from\n    " << meshDict.filePath() << endl;

  blockMesh blocks(meshDict, regionName);

  Info << nl << "Creating polyMesh from blockMesh" << endl;

  word defaultFacesName = "defaultFaces";
  word defaultFacesType = emptyPolyPatch::typeName;

  polyMesh mesh(IOobject(regionName, runTime->constant(), *runTime),
                Foam::pointField(blocks.points()), blocks.cells(),
                blocks.patches(), blocks.patchNames(), blocks.patchDicts(),
                defaultFacesName, defaultFacesType);

  // Disabling mergePairs feature for now
  /*// Read in a list of dictionaries for the merge patch pairs
  if (meshDict.found("mergePatchPairs"))
  {
      List<Pair<word>> mergePatchPairs
      (
          meshDict.lookup("mergePatchPairs")
      );

      #include "mergePatchPairs.H"
  }
  else
  {
      Info<< nl << "There are no merge patch pairs edges" << endl;
  }*/

  label nZones = blocks.numZonedBlocks();

  if (nZones > 0) {
    Info << nl << "Adding cell zones" << endl;

    // Map from zoneName to cellZone index
    HashTable<label> zoneMap(nZones);

    // Cells per zone.
    List<DynamicList<label>> zoneCells(nZones);

    // Running cell counter
    label celli = 0;

    // Largest zone so far
    label freeZoneI = 0;

    forAll(blocks, blockI) {
      const block& b = blocks[blockI];

      const List<FixedList<label, 8>> blockCells = b.cells();

      const word& zoneName = b.zoneName();

      if (zoneName.size()) {
        HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

        label zoneI;

        if (iter == zoneMap.end()) {
          zoneI = freeZoneI++;

          Info << "    " << zoneI << '\t' << zoneName << endl;

          zoneMap.insert(zoneName, zoneI);
        } else {
          zoneI = iter();
        }

        forAll(blockCells, i) { zoneCells[zoneI].append(celli++); }

      } else {
        celli += b.cells().size();
      }
    }

    List<cellZone*> cz(zoneMap.size());

    Info << nl << "Writing cell zones as cellSets" << endl;

    forAllConstIter(HashTable<label>, zoneMap, iter) {
      label zoneI = iter();

      cz[zoneI] = new cellZone(iter.key(), zoneCells[zoneI].shrink(), zoneI,
                               mesh.cellZones());

      // Write as cellSet for ease of processing
      cellSet cset(mesh, iter.key(), zoneCells[zoneI].shrink());
      cset.write();
    }

    mesh.pointZones().setSize(0);
    mesh.faceZones().setSize(0);
    mesh.cellZones().setSize(0);
    mesh.addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
  }

  // Set the precision of the points data to 10
  IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

  if (write) {
    Info << nl << "Writing polyMesh" << endl;
    mesh.removeFiles();
    if (!mesh.write()) {
      FatalErrorInFunction << "Failed writing polyMesh." << exit(FatalError);
    }
  }

  //
  // writes some information
  //
  {
    const polyPatchList& patches = mesh.boundaryMesh();

    Info << "----------------" << nl << "Mesh Information" << nl
         << "----------------" << nl << "  "
         << "boundingBox: " << boundBox(mesh.points()) << nl << "  "
         << "nPoints: " << mesh.nPoints() << nl << "  "
         << "nCells: " << mesh.nCells() << nl << "  "
         << "nFaces: " << mesh.nFaces() << nl << "  "
         << "nInternalFaces: " << mesh.nInternalFaces() << nl;

    Info << "----------------" << nl << "Patches" << nl << "----------------"
         << nl;

    forAll(patches, patchi) {
      const polyPatch& p = patches[patchi];

      Info << "  "
           << "patch " << patchi << " (start: " << p.start()
           << " size: " << p.size() << ") name: " << p.name() << nl;
    }
  }

  Info << "\nEnd\n" << endl;

  // Copy polyMesh into fvMesh now and return that
  Foam::PtrList<Foam::polyPatch> patchPtrList(0);
  const polyBoundaryMesh& patches2 = mesh.boundaryMesh();

  for (int i = 0; i < patches2.size(); i++) {
    Foam::polyPatch* tmpPtch = new Foam::polyPatch(patches2[i]);
    patchPtrList.append(tmpPtch);
  }

  // Get pointField (points), faceList (faces), and cellList (cells).
  const Foam::cellList cellLst = mesh.cells();
  const Foam::faceList faceLst = mesh.faces();
  const Foam::pointField ptsAll = mesh.points();

  Foam::pointField pointData2(mesh.nPoints());
  Foam::cellList cellLst2(0);
  cellLst2 = cellLst;
  Foam::faceList faceLst2(0);
  faceLst2 = faceLst;
  pointData2 = ptsAll;

  auto* fm = new Foam::fvMesh(
      IOobject(regionName, runTime->timeName(), *runTime,
               Foam::IOobject::NO_READ, Foam::IOobject::AUTO_WRITE),
      std::move(pointData2), std::move(faceLst2), std::move(cellLst2));

  fm->addFvPatches(patchPtrList, true);
  return fm;
}