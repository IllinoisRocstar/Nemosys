#ifndef PACKMESHDRIVER_H
#define PACKMESHDRIVER_H

#include "NemDriver.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

#ifdef HAVE_CFMSH
// Friend classes for parameters
class MeshManipulationFoamParams;
class cfmeshParams;
class blockMeshParams;
class snappymeshParams;
class MeshManipulationFoam;

#endif

/**
  @brief This class drives the flow of pack meshing by using parameters provided
  by users. It supports two different meshing pipelines. First being hexahedral
  conformal meshing and second is tetrahedral periodic as well as non-priodic
  meshing.

  Hexahedral meshing pipeline creates different objects for different meshing
  engines in NEMoSys, as well as objects of mesh manipulation utilities to
  complete the meshing workflow.
  ## Pack Meshing Workflow for Hexahedral Mesh
  > * Geometry generation using RocPack output file.
  > * Generation of STL file for workflow downstream.
  > * Generation of background mesh for snappyHexMesh
  > * Snapping of pack surfaces onto background volume mesh to create conformal
interfaces. > * Splitting multiple pack regions into different cellZones
(volumes). > * Merging all different volumes into two different foam mesh (1.
Packs, 2. Surrounding) > * Creating patch for Packs and Surrounding regions. > *
Conversion of foam mesh to VTK. > * Provides mesh quality statistics to users in
a text file.

  Tetrahedral periodic meshing pipeline operates within rocPack class and uses
  Netgen mesh optimization for final mesh quality improvement if needed.
  ## Pack Meshing Workflow for Tetrahedral Mesh
  > * Geometry generation using RocPack output file.
  > * Enforcing periodicity in geometry for periodic mesh
  > * Mapping periodic surfaces and points with their counterparts.
  > * Meshing the geometry tetrahedrally while enforcing periodic boundary
  boundary conditions.
  > * Exporting final mesh with or without optimization.

**/
class PackMeshDriver : public NemDriver {
 public:
#ifdef HAVE_CFMSH
  /**
      @brief PackMeshDriver constructor for hexahedral meshing
      @param ifname Input surface file name
      @param _mparams MeshManipulationFoamParams object
      @param _cfparams cfmeshParams object
      @param _snappyparams snappymeshParams object
      @param _bmparams blockMeshParams object
      @param ofname1 Pack mesh output filename.extension
      @param useRocpack Boolean for using rocpack to generate an STL
      @param ofname2 Surrounding mesh output filename.extension
  **/
  PackMeshDriver(const std::string &ifname,
                 MeshManipulationFoamParams *_mparams, cfmeshParams *_cfparams,
                 snappymeshParams *_snappyparams, blockMeshParams *_bmparams,
                 const std::string &ofname1, const std::string &ofname2,
                 const bool &useRocpack);
#endif

  /**
      @brief PackMeshDriver alternate constructor for tetrahedral meshing
      @param ifname Input file name (Output from Rocpack)
      @param ofname Output file name (.stl, .vtu, .msh are supported)
      @param scaleVolumes Boolean for scaling the volumes
      @param scaleValue Value of scaling between 0 to 1 with 1 being no scaling.
      @param meshSize Mesh size
      @param rmvBndryPacks Boolean for removal of volumes intersecting boundary
      @param setPeriodicGeom Boolean for setting periodic geometry
      @param setPeriodicMesh Boolean for setting periodic mesh
      @param enable2PhysGrps Boolean for adding two physical groups
      @param enableMultiPhysGrps Bool for adding several physical groups
      @param wantGeometryOnly Bool for getting stl file only
      @param createCohesive Boolean for enabling cohesive elements creation
      @param enablePatches Boolean for enabling patches.
      @param transferMesh Parameters for translating mesh from origin.
      @param customDomain Boolean for users to specify custom domain size.
      @param domainBounds User-defined domain size.
      @param mshAlgorithm Meshing algorithm of choice.
      @param enableDefaultOutput Enable default output formats
      @param enablePhysGrpPerShape Enables physical group per shape
      @param refineLevel Refinement applied to original mesh
      @param upperThreshold Upper threshold for filtering (w.r.t mean)
      @param lowerThreshold Lower threshold for filtering (w.r.t mean)
      @param preserveSize Boolean for preserving pack sizes instead of
             packing fraction
  **/
  PackMeshDriver(const std::string &ifname, const std::string &ofname,
                 const bool &scaleVolumes, const double &scaleValue,
                 const double &meshSize, const bool &rmvBndryPacks,
                 const bool &setPeriodicGeom, const bool &setPeriodicMesh,
                 const bool &enable2PhysGrps, const bool &enableMultiPhysGrps,
                 const bool &wantGeometryOnly, const bool &createCohesive,
                 const bool &enablePatches,
                 const std::vector<double> &transferMesh,
                 const bool &customDomain,
                 const std::vector<double> &domainBounds,
                 const int &mshAlgorithm, const bool &enableDefaultOutput,
                 const bool &enablePhysGrpPerShape, const int &refineLevel,
                 const double &upperThreshold, const double &lowerThreshold,
                 const bool &preserveSize);

  /** @brief PackMeshDriver standard destructor
   **/
  ~PackMeshDriver();

  // static factory method
  /** @brief reads JSON object
      @param inputjson JSON object created using user provided input file
  **/
  static PackMeshDriver *readJSON(const jsoncons::json inputjson);

#ifdef HAVE_CFMSH
  /** @brief Input parsing module for hex mesh. Assigned user defined parameters
             in input file to variables.
      @param ifname Input surface file name
      @param ofname1 Pack mesh output filename.extension
      @param ofname2 Surrounding mesh output filename.extension
      @param useRocpack Boolean for using rocpack to generate an STL
      @param inputjson JSON object created using user provided input file
  **/
  static PackMeshDriver *readJSON(const std::string &ifname,
                                  const std::string &ofname1,
                                  const std::string &ofname2,
                                  const bool useRocpack,
                                  const jsoncons::json inputjson);
#endif

  /** @brief Input parsing module for tetrahedral periodic/non-periodic mesh
      @param ifname Input file name
      @param ofname Output file name
      @param inputjson JSON object created using user provided input file
  **/
  static PackMeshDriver *readJSON(const std::string &ifname,
                                  const std::string &ofname,
                                  const jsoncons::json inputjson);

 private:
};

#endif
