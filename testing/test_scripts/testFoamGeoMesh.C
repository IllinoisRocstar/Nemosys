#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#  define _USE_MATH_DEFINES
#endif
#include <gtest/gtest.h>

#include <Mesh/foamGeoMesh.H>

#include <argList.H>
#include <cellModeller.H>
#include <fileName.H>
#include <fvMesh.H>

std::string fileNameHexRead;
std::string fileNameTetRead;
std::string fileNameHexWrite;
std::string fileNameTetWrite;
Foam::fvMesh *fm;

Foam::fvMesh* getFoamMesh() {
  Foam::dictionary *ControlDict1_ = new Foam::dictionary("controlDict");
  Foam::dictionary fmfle("FoamFile");
  fmfle.add("version", "2.0");
  fmfle.add("format", "ascii");
  fmfle.add("class", "dictionary");
  fmfle.add("location", "\"system\"");
  fmfle.add("object", "controlDict");
  ControlDict1_->add("FoamFile",fmfle);
  ControlDict1_->add("deltaT",1);
  ControlDict1_->add("startTime",0);
  ControlDict1_->add("writeInterval",1);

  // Create time class without reading controlDict
  Foam::Time *runTime_ = new Foam::Time(ControlDict1_, ".", ".");
  Foam::argList::noParallel();

  Foam::pointField pointData(0);
  Foam::cellShapeList cellShapeData(0);

  auto tmpfm = std::unique_ptr<Foam::polyMesh>(new Foam::polyMesh(
    Foam::IOobject
    (
      Foam::polyMesh::defaultRegion,
      runTime_->constant(),
      *runTime_
    ),
    std::move(pointData),   // Vertices
    cellShapeData,               // Cell shape and points
    Foam::faceListList(),        // Boundary faces
    Foam::wordList(),            // Boundary Patch Names
    Foam::PtrList<Foam::dictionary>(), // Boundary Dicts
    "defaultPatch",              // Default Patch Name
    Foam::polyPatch::typeName    // Default Patch Type
  ));

  Foam::fvMesh* fm = dynamic_cast<Foam::fvMesh*>(tmpfm.get());
  return fm;
}

bool comparePatches(Foam::fileName &a, Foam::fileName &b) {
  std::vector<std::string> onePatches;
  std::vector<int> onePts;
  std::vector<std::string> twoPatches;
  std::vector<int> twoPts;

  Foam::dictionary *ControlDict1_ = new Foam::dictionary("controlDict");
  Foam::dictionary fmfle("FoamFile");
  fmfle.add("version", "2.0");
  fmfle.add("format", "ascii");
  fmfle.add("class", "dictionary");
  fmfle.add("location", "\"system\"");
  fmfle.add("object", "controlDict");
  ControlDict1_->add("FoamFile",fmfle);
  ControlDict1_->add("deltaT",1);
  ControlDict1_->add("startTime",0);
  ControlDict1_->add("writeInterval",1);

  // Create time class without reading controlDict
  Foam::Time *runTime_ = new Foam::Time(ControlDict1_, ".", ".");
  Foam::argList::noParallel();

  // a
  {
    auto *msh = new Foam::fvMesh(
      Foam::IOobject
      (
          a,
          runTime_->constant(),
          *runTime_,
          Foam::IOobject::MUST_READ,
          Foam::IOobject::AUTO_WRITE
      )
    );

    const Foam::polyBoundaryMesh& patches = msh->boundaryMesh();
    Foam::wordList patchNames = patches.names();
    for (int i=0; i<patchNames.size(); i++) onePatches.push_back(patchNames[i]);
    for (int i=0; i<patches.size(); i++) onePts.push_back(patches[i].nPoints());
  }

  // b
  {
    auto *msh = new Foam::fvMesh(
      Foam::IOobject
      (
          b,
          runTime_->constant(),
          *runTime_,
          Foam::IOobject::MUST_READ,
          Foam::IOobject::AUTO_WRITE
      )
    );

    const Foam::polyBoundaryMesh& patches = msh->boundaryMesh();
    Foam::wordList patchNames = patches.names();
    for (int i=0; i<patchNames.size(); i++) twoPatches.push_back(patchNames[i]);
    for (int i=0; i<patches.size(); i++) twoPts.push_back(patches[i].nPoints());
  }

  // Compare two arrays
  for (int i=0; i<onePatches.size(); i++)
    if (onePatches[i] != twoPatches[i])
      return false;
  for (int i=0; i<onePts.size(); i++)
    if (onePts[i] != twoPts[i])
      return false;

  return true;
}

TEST(foamGeoMesh, ConstructorDefault) { NEM::MSH::foamGeoMesh fgm{}; }

TEST(foamGeoMesh, ConstructorFvMesh) {
  fm = getFoamMesh();
  NEM::MSH::foamGeoMesh fgm{fm};
}

TEST(foamGeoMesh, ConstructorFvMeshWithGeo) {
  fm = getFoamMesh();
  NEM::MSH::foamGeoMesh fgm{fm, "GeoIds"};
}

TEST(foamGeoMesh, ReadWriteHexMesh) {
  NEM::MSH::foamGeoMesh *fgm = NEM::MSH::foamGeoMesh::Read(fileNameHexRead,"PhysIds");
  fgm->write(fileNameHexWrite);
  fgm->Delete();
}

TEST(foamGeoMesh, ReadWriteTetMesh) {
  NEM::MSH::foamGeoMesh *fgm = NEM::MSH::foamGeoMesh::Read(fileNameTetRead,"PhysIds");
  fgm->write(fileNameTetWrite);
  fgm->Delete();
}

TEST(foamGeoMesh, GM2FoamTest) {
  NEM::MSH::foamGeoMesh *fgm = NEM::MSH::foamGeoMesh::Read(fileNameHexRead,"PhysIds");
  fgm->resetNative();
  fgm->write("reset_mesh");
  fgm->Delete();
}

TEST(foamGeoMesh, checkPatches) {
  Foam::fileName refFile = fileNameHexRead;
  Foam::fileName testFile = "reset_mesh";
  bool tstBool = comparePatches(refFile, testFile);

  refFile = fileNameHexRead;
  testFile = fileNameHexWrite;
  bool tstBool2 = comparePatches(refFile, testFile);

  EXPECT_EQ(tstBool, 1);
  EXPECT_EQ(tstBool2, 1);
}


int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  assert(argc == 5);
  fileNameHexRead  = argv[1];
  fileNameTetRead  = argv[2];
  fileNameHexWrite = argv[3];
  fileNameTetWrite = argv[4];

  int res = RUN_ALL_TESTS();

  return res;
}
