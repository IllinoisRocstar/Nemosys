#ifndef PACKMESHDRIVER_H
#define PACKMESHDRIVER_H

#include "NemDriver.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

// Friend classes for parameters
class MeshManipulationFoamParams;
class cfmeshParams;
class blockMeshParams;
class snappymeshParams;
class MeshManipulationFoam;


/**
  @brief This class drives the flow of pack meshing by using parameters provided
  by users. It creates different objects for different meshing engines in NEMoSys,
  as well as objects of mesh manipulation utilities to complete the meshing workflow.
  ## Pack Meshing Workflow
  > * Hexahedral mesh generation of pack geometry using CfMesh.
  > * Extraction of pack mesh surfaces.
  > * Generation of background mesh for snappyHexMesh
  > * Snapping of pack surfaces onto background volume mesh to create conformal interfaces.
  > * Splitting multiple pack regions into different cellZones (volumes).
  > * Merging all different volumes into two different foam mesh (1. Packs, 2. Surrounding)
  > * Creating patch for Packs and Surrounding regions.
  > * Conversion of foam mesh to VTK.
  > * Provides mesh quality statistics to users in a text file.
**/
class PackMeshDriver : public NemDriver
{
  public:
  /** 
      @brief PackMeshDriver constructor
      @param ifname Input surface file name
      @param _mparams MeshManipulationFoamParams object
      @param _cfparams cfmeshParams object
      @param _snappyparams snappymeshParams object
      @param _bmparams blockMeshParams object
      @param ofname1 Pack mesh output filename.extension
      @param ofname2 Surrounding mesh output filename.extension
  **/
  PackMeshDriver(const std::string& ifname, 
                MeshManipulationFoamParams* _mparams,
                cfmeshParams* _cfparams, snappymeshParams* _snappyparams,
                blockMeshParams* _bmparams, const std::string& ofname1,
                const std::string& ofname2);

  /** @brief PackMeshDriver standard destructor
  **/
  ~PackMeshDriver();
    
  // static factory method
  /** @brief reads JSON object
      @param inputjson JSON object created using user provided input file
  **/
  static PackMeshDriver* readJSON(const jsoncons::json inputjson);

  /** @brief Input parsing module. Assigned user defined parameters in input 
            file to variables.
      @param ifname Input surface file name
      @param ofname1 Pack mesh output filename.extension
      @param ofname2 Surrounding mesh output filename.extension
      @param inputjson JSON object created using user provided input file
  **/
  static PackMeshDriver* readJSON(const std::string& ifname,
                                  const std::string& ofname1,
                                  const std::string& ofname2,
                                  const jsoncons::json inputjson);

  private:
  // Nothing for now

};

#endif
