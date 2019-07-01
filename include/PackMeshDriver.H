#ifndef PACKMESHDRIVER_H
#define PACKMESHDRIVER_H

#include "NemDriver.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
/*
  This class drives the flow of pack meshing by using parameters provided
  by users. It creates different objects for different meshing engines,
  as well as object of mesh manipulation utility to faciilitate the meshing.
*/

// Friend classes for parameters
class MeshManipulationFoamParams;
class cfmeshParams;
class blockMeshParams;
class snappymeshParams;
class MeshManipulationFoam;

class PackMeshDriver : public NemDriver
{
  public:
  // standard constructor
  PackMeshDriver(const std::string& ifname, 
                MeshManipulationFoamParams* _mparams,
                cfmeshParams* _cfparams, snappymeshParams* _snappyparams,
                blockMeshParams* _bmparams, const std::string& ofname1,
                const std::string& ofname2);

  // default descrtuctor
  ~PackMeshDriver();
    
  // static factory method
  static PackMeshDriver* readJSON(const jsoncons::json inputjson);
  static PackMeshDriver* readJSON(const std::string& ifname,
                                  const std::string& ofname1,
                                  const std::string& ofname2,
                                  const jsoncons::json inputjson);

  private:
  // Nothing for now

};

#endif
