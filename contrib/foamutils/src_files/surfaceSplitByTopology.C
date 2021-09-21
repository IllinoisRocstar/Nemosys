#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "getDicts.H"
#include "surfaceSplitByTopology.H"

using namespace Foam;

surfaceSplitByTopology::surfaceSplitByTopology() {}

surfaceSplitByTopology::~surfaceSplitByTopology() {}

int surfaceSplitByTopology::execute(const Foam::fileName surfFileName,
                                    const Foam::fileName outFileName) {
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);
  auto fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  auto fvSolution_ = initFoam->createFvSolution(writeDicts);

  // Create time class without reading controlDict
  Foam::Time runTime(controlDict_.get(), ".", ".");
  Foam::argList::noParallel();

  fileName outFileBaseName = outFileName.lessExt();
  word outExtension = outFileName.ext();

  // Load surface
  triSurface surf(surfFileName);

  bool anyZoneRemoved = false;

  label iterationNo = 0;
  label iterationLimit = 10;

  Info << "Splitting off baffle parts " << endl;

  do {
    anyZoneRemoved = false;

    labelList faceZone;

    const labelListList& edFaces = surf.edgeFaces();
    const labelListList& faceEds = surf.faceEdges();

    boolList multipleEdges(edFaces.size(), false);

    forAll(multipleEdges, i) {
      if (edFaces[i].size() > 2) {
        multipleEdges[i] = true;
      }
    }

    label nZones = surf.markZones(multipleEdges, faceZone);

    if (nZones < 2) {
      break;
    }

    boolList nonBaffle(faceZone.size(), true);
    boolList baffle(faceZone.size(), true);
    labelList pointMap;
    labelList faceMap;

    for (label z = 0; z < nZones; z++) {
      bool keepZone = true;

      forAll(faceZone, f) {
        if (faceZone[f] == z) {
          forAll(faceEds[f], fe) {
            if (edFaces[faceEds[f][fe]].size() < 2) {
              keepZone = false;

              anyZoneRemoved = true;

              break;
            }
          }
        }

        if (!keepZone) {
          break;
        }
      }

      forAll(faceZone, f) {
        if (faceZone[f] == z) {
          nonBaffle[f] = keepZone;
          baffle[f] = !keepZone;
        }
      }
    }

    Info << "    Iteration " << iterationNo << endl;

    triSurface baffleSurf = surf.subsetMesh(baffle, pointMap, faceMap);

    if (baffleSurf.size()) {
      fileName bafflePartFileName = outFileBaseName + "_bafflePart_" +
                                    name(iterationNo) + "." + outExtension;

      Info << "    Writing baffle part to " << bafflePartFileName << endl;

      baffleSurf.write(bafflePartFileName);
    }

    surf = surf.subsetMesh(nonBaffle, pointMap, faceMap);

    if (iterationNo == iterationLimit) {
      WarningInFunction << "Iteration limit of " << iterationLimit << "reached"
                        << endl;
    }

    iterationNo++;

  } while (anyZoneRemoved && iterationNo < iterationLimit);

  Info << "Writing new surface to " << outFileName << endl;

  surf.write(outFileName);

  labelList faceZone;

  const labelListList& edFaces = surf.edgeFaces();

  boolList multipleEdges(edFaces.size(), false);

  forAll(multipleEdges, i) {
    if (edFaces[i].size() > 2) {
      multipleEdges[i] = true;
    }
  }

  label nZones = surf.markZones(multipleEdges, faceZone);

  int nofPacks = nZones;

  // nZones is number of pack domains. it will be great to get that out for
  // mergeMeshes. Also, surface renumbering is needed as output STL

  /*Info<< "Splitting remaining multiply connected parts" << endl;

  for (label z = 0; z < nZones; z++)
  {

      boolList include(faceZone.size(), false);
      labelList pointMap;
      labelList faceMap;

      forAll(faceZone, f)
      {
          if (faceZone[f] == z)
          {
              include[f] = true;
          }
      }

      triSurface zoneSurf = surf.subsetMesh(include, pointMap, faceMap);


      fileName remainingPartFileName =
          outFileBaseName
        + "_multiplePart_"
        + name(z)
        + "." + outExtension;

      Info<< "    Writing mulitple part "
          << z << " to " << remainingPartFileName << endl;

      zoneSurf.write(remainingPartFileName);
  }*/

  Info << "End\n" << endl;

  return nofPacks;
}