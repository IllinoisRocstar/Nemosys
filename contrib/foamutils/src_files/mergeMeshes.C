#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "getDicts.H"
#include "mergeMeshes.H"

using namespace Foam;

mergeMeshes::mergeMeshes() {}

mergeMeshes::~mergeMeshes() {}

void mergeMeshes::execute(const std::string& mainCase,
                          const std::string& addCasePath,
                          const std::vector<std::string>& addCases,
                          const int& nDomains, const int& dirStat) {
  bool writeDicts = false;
  std::unique_ptr<getDicts> initFoam;
  initFoam = std::unique_ptr<getDicts>(new getDicts());
  auto controlDict_ = initFoam->createControlDict(writeDicts);
  auto fvSchemes_ = initFoam->createFvSchemes(writeDicts);
  auto fvSolution_ = initFoam->createFvSolution(writeDicts);

  // Create time class without reading controlDict
  Foam::Time runTime(controlDict_.get(), ".", ".");
  Foam::argList::noParallel();

  // Main for loop. It loops through all the slave regions and adds them to
  // master region one by one. At the end, master region directory will have
  // all the slave regions added. Keep overwrite boolean true for this.
  for (int j = 0; j < (nDomains - 1); j++) {
    const bool overwrite = true;
    word masterRegion = polyMesh::defaultRegion;

    fileName masterCase = mainCase;
    if (nDomains == 2) {
      if (dirStat == 1)
        masterRegion = "domain2";
      else
        masterRegion = "domain1";
    } else
      masterRegion = "domain1";

    fileName addCase = addCasePath;
    word addRegion = polyMesh::defaultRegion;
    addRegion = addCases[j];

    getRootCase(masterCase);
    getRootCase(addCase);

    Info << "Master:      " << masterCase << "  region " << masterRegion << nl
         << "mesh to add: " << addCase << "  region " << addRegion << endl;

    const fileName masterCasePath = masterCase.path();
    const fileName masterCaseName = masterCase.name();

    Time runTimeMaster(Time::controlDictName, masterCasePath, masterCaseName);
    runTimeMaster.functionObjects().off();

    const fileName addCasePath = addCase.path();
    const fileName addCaseName = addCase.name();

    Time runTimeToAdd(Time::controlDictName, addCasePath, addCaseName);
    runTimeToAdd.functionObjects().off();

    Info << "Reading master mesh for time = " << runTimeMaster.timeName() << nl;

    Info << "Create mesh\n" << endl;
    Foam::mergePolyMesh masterMesh(
        IOobject(masterRegion, runTimeMaster.timeName(), runTimeMaster));
    const word oldInstance = masterMesh.pointsInstance();

    Info << "Reading mesh to add for time = " << runTimeToAdd.timeName() << nl;

    Info << "Create mesh\n" << endl;
    polyMesh meshToAdd(
        IOobject(addRegion, runTimeToAdd.timeName(), runTimeToAdd));

    if (!overwrite) {
      runTimeMaster++;
    }

    Info << "Writing combined mesh to " << runTimeMaster.timeName() << endl;

    masterMesh.addMesh(meshToAdd);
    masterMesh.merge();

    if (overwrite) {
      masterMesh.setInstance(oldInstance);
    }

    masterMesh.write();
  }

  Info << "End\n" << endl;
}

void mergeMeshes::getRootCase(Foam::fileName& casePath) {
  using namespace Foam;
  casePath.clean();

  if (casePath.empty() || casePath == ".") {
    // handle degenerate form and '.'
    casePath = cwd();
  } else if (casePath[0] != '/' && casePath.name() == "..") {
    // avoid relative cases ending in '..' - makes for very ugly names
    casePath = cwd() / casePath;
    casePath.clean();
  }
}