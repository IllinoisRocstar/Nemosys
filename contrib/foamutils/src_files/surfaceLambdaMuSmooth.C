#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "getDicts.H"
#include "surfaceLambdaMuSmooth.H"

using namespace Foam;

surfaceLambdaMuSmooth::surfaceLambdaMuSmooth() {}

surfaceLambdaMuSmooth::~surfaceLambdaMuSmooth() {}

void surfaceLambdaMuSmooth::execute(const Foam::fileName surfFileName,
                                    const Foam::fileName outFileName,
                                    const Foam::scalar lambda,
                                    const Foam::scalar mu,
                                    const Foam::label iters,
                                    const bool addFtrFl) {
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);
  auto fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  auto fvSolution_ = initFoam->createFvSolution(writeDicts);

  // Create time class without reading controlDict
  Foam::Time runTime(controlDict_.get(), ".", ".");
  Foam::argList::noParallel();

  if (lambda < 0 || lambda > 1) {
    FatalErrorInFunction
        << lambda << endl
        << "0: no change   1: move vertices to average of neighbours"
        << exit(FatalError);
  }
  if (mu < 0 || mu > 1) {
    FatalErrorInFunction
        << mu << endl
        << "0: no change   1: move vertices to average of neighbours"
        << exit(FatalError);
  }
  Info << "lambda      : " << lambda << nl << "mu          : " << mu << nl
       << "Iters       : " << iters << nl << "Reading surface from "
       << surfFileName << " ..." << endl;

  meshedSurface surf1(surfFileName);

  Info << "Faces       : " << surf1.size() << nl
       << "Vertices    : " << surf1.nPoints() << nl
       << "Bounding Box: " << boundBox(surf1.localPoints()) << endl;

  bitSet fixedPoints(surf1.localPoints().size(), false);

  if (addFtrFl) {
    const fileName featureFileName("ftrEdge.ftr");
    Info << "Reading features from " << featureFileName << " ..." << endl;

    edgeMesh feMesh(featureFileName);

    getFixedPoints(feMesh, surf1.localPoints(), fixedPoints);

    Info << "Number of fixed points on surface = " << fixedPoints.count()
         << endl;
  }

  pointField newPoints(surf1.localPoints());

  for (label iter = 0; iter < iters; iter++) {
    // Lambda
    {
      pointField newLocalPoints((1 - lambda) * surf1.localPoints() +
                                lambda * avg(surf1, fixedPoints));

      pointField newPoints(surf1.points());
      UIndirectList<point>(newPoints, surf1.meshPoints()) = newLocalPoints;

      surf1.movePoints(newPoints);
    }

    // Mu
    if (mu != 0) {
      pointField newLocalPoints((1 + mu) * surf1.localPoints() -
                                mu * avg(surf1, fixedPoints));

      pointField newPoints(surf1.points());
      UIndirectList<point>(newPoints, surf1.meshPoints()) = newLocalPoints;

      surf1.movePoints(newPoints);
    }
  }

  Info << "Writing surface to " << outFileName << " ..." << endl;
  surf1.write(outFileName);

  Info << "End\n" << endl;
}

Foam::tmp<Foam::pointField> surfaceLambdaMuSmooth::avg(
    const Foam::meshedSurface& s, const Foam::bitSet& fixedPoints) {
  const labelListList& pointEdges = s.pointEdges();

  tmp<pointField> tavg(new pointField(s.nPoints(), Zero));
  pointField& avg = tavg.ref();

  forAll(pointEdges, vertI) {
    vector& avgPos = avg[vertI];

    if (fixedPoints[vertI]) {
      avgPos = s.localPoints()[vertI];
    } else {
      const labelList& pEdges = pointEdges[vertI];

      forAll(pEdges, myEdgeI) {
        const edge& e = s.edges()[pEdges[myEdgeI]];

        label otherVertI = e.otherVertex(vertI);

        avgPos += s.localPoints()[otherVertI];
      }

      avgPos /= pEdges.size();
    }
  }

  return tavg;
}

// SurfaceLambdaMuSmooth
void surfaceLambdaMuSmooth::getFixedPoints(const Foam::edgeMesh& feMesh,
                                           const Foam::pointField& points,
                                           Foam::bitSet& fixedPoints) {
  scalarList matchDistance(feMesh.points().size(), 1e-1);
  labelList from0To1;

  bool matchedAll =
      matchPoints(feMesh.points(), points, matchDistance, false, from0To1);

  if (!matchedAll) {
    WarningInFunction
        << "Did not match all feature points to points on the surface" << endl;
  }

  forAll(from0To1, fpI) {
    if (from0To1[fpI] != -1) {
      fixedPoints[from0To1[fpI]] = true;
    }
  }
}