#ifndef NEMOSYS_CFMESHGEN_H_
#define NEMOSYS_CFMESHGEN_H_

#ifdef HAVE_CFMSH

#include "MeshGeneration/meshGen.H"

// OpenFOAM headers
#include <argList.H>
#include <fvOptions.H>

class cfmeshParams;

class cfmeshGen : public meshGen {
 public:
  cfmeshGen();
  explicit cfmeshGen(cfmeshParams *params);
  ~cfmeshGen() override;

  // cfmesh mesh creation
 public:
  // create a 3D polygonal discretization from STL file
  int createMeshFromSTL(const char *fname) override;

  // internal management
 private:
  // pre-processor
  void initialize();

  // create dictionaries needed for openFoam
  void createMshDict(const bool &write);

  // surface edge extraction
  int surfaceFeatureEdgeDetect();

  // mesh quality improvement
  int improveMeshQuality();

  // class data and management flags
 private:
  bool defaults;
  std::string caseName{"_cfmesh"};  // internal information
  cfmeshParams *params_;

  // openfoam data structure
 private:
  std::unique_ptr<Foam::dictionary> controlDict_;
  std::unique_ptr<Foam::dictionary> fvSchemes_;
  std::unique_ptr<Foam::dictionary> fvSolution_;
  std::unique_ptr<Foam::Time> runTime_;
  std::unique_ptr<Foam::fvMesh> fmesh_;
  std::unique_ptr<Foam::dictionary> meshDict_;
};

#endif

#endif  // NEMOSYS_CFMESHGEN_H_
