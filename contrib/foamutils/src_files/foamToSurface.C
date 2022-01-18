#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "foamToSurface.H"
#include "getDicts.H"

using namespace Foam;

foamToSurface::foamToSurface() {}

foamToSurface::~foamToSurface() {}

void foamToSurface::execute(const std::string &outputfile) {
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);
  auto fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  auto fvSolution_ = initFoam->createFvSolution(writeDicts);

  // Create time class without reading controlDict
  Foam::Time runTime(controlDict_.get(), ".", ".");
  Foam::argList::noParallel();

  fileName exportName = outputfile;

  scalar scaleFactor = 0;
  const bool doTriangulate = true;

  fileName exportBase = exportName.lessExt();
  word exportExt = exportName.ext();

  polyMesh mesh(IOobject(Foam::polyMesh::defaultRegion, runTime.timeName(),
                         runTime, Foam::IOobject::MUST_READ));

  // polyMesh::readUpdateState state = mesh.readUpdate();
  mesh.readUpdate();
  exportName = exportBase + "." + exportExt;

  meshedSurface surf(mesh.boundaryMesh());
  surf.scalePoints(scaleFactor);

  Info << "writing " << exportName;

  if (doTriangulate) {
    Info << " triangulated";
    surf.triangulate();
  }

  if (scaleFactor <= 0) {
    Info << " without scaling" << endl;
  } else {
    Info << " with scaling " << scaleFactor << endl;
  }

  // Add an error handler for directory check before writing
  surf.write(exportName);

  Info << nl << endl;

  Info << "End\n" << endl;
}