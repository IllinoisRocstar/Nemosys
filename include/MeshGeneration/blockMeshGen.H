#ifndef NEMOSYS_BLOCKMESHGEN_H_
#define NEMOSYS_BLOCKMESHGEN_H_

#ifdef HAVE_CFMSH

// Nemosys Headers
#include "MeshGeneration/meshGen.H"
#include "Mesh/geoMeshBase.H"
#include "Mesh/foamGeoMesh.H"

// Foam headers
#include <argList.H>
#include <fvOptions.H>

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

class blockMeshParams;  // Friend Class

/**
  @brief blockMeshGen <-- meshGen <-- meshBase
        This class incorporates mesh generating method of blockMesh utility.
        Currently, we are supporting automatic full-hex mesh generation
        Box, Sphere, Cylinder, and Tapered Cone geometries. Output mesh is
        written in OpenFOAM polyMesh format.
**/
class blockMeshGen : public meshGen {
  // --- constructors and destructor
 public:
  /** @brief blockMeshGen standard constructor
   **/
  blockMeshGen();  // Default constructor

  /** @brief blockMeshGen alternate constructor. Uses user-defined parameters
            to perform requested meshing operation.
      @param params blockMeshParams object
  **/
  blockMeshGen(blockMeshParams *params);  // Constructor with parameters

  /** @brief blockMeshGen standard desctructor
   **/
  ~blockMeshGen();  // Class destructor

  // -- meshGen methods override
 public:
  /** @brief Generates mesh and returns VTK database
      @param fname Input surface file name
      @return meshBase mesh dataset
  **/
  int createMeshFromSTL(const char *fname);

  // -- Internal
 private:
  /**
      @brief Initializes OpenFOAM args and runtime
  **/
  void initialize();

  /**
      @brief Creates blockMeshDict from user arguments
  **/
  void createBlockMshDict(const bool &write);

  // --- Internal Use
 private:
  bool defaults;  // Default parameters boolean

  /** @brief Definition of parameters pointer
   **/
  blockMeshParams *params_;  // Definition of parameters pointer

  // --- OpenFOAM Data Structure
 private:
  std::unique_ptr<Foam::Time> runTime_;
  std::unique_ptr<Foam::fvMesh> fmesh_;
  std::unique_ptr<Foam::dictionary> blockMshDict_;
  std::unique_ptr<Foam::dictionary> controlDict_;
  std::unique_ptr<Foam::dictionary> fvSchemes_;
  std::unique_ptr<Foam::dictionary> fvSolution_;
};

#endif

#endif  // NEMOSYS_BLOCKMESHGEN_H_
