#ifndef NEMDRIVERS_H
#define NEMDRIVERS_H

#include <meshUser.H>
#include <jsoncons/json.hpp>


/* ---------------------------------------------------- 
TODO: This class will contain drivers for several utils:
    1) Mesh Generation :
      - ng_meshing_params
      - stl file
    2) Data Transfer:
      - mesh file with data (source)
      - mesh file without data (target)
      - name of array to transfer
        OR no name (transfers all)
      - method to use for transfer (currently only "FE")
    3) Adaptive Mesh Refinement:
      - mesh file with data
      - method of refinement ("gradient", "value", "error")
      - name of array to use for size field gen (w/ method)
    4) Uniform Mesh Refinement:
      - mesh file w/ or w/o data
      - scaling factor for edges   
    5) Format conversion:
      - vt* to .msh, .vol etc.
      - .msh to vtu
      - .vol to vtu 

   -------------------------------------------------- */


using jsoncons::json;

class NemDriver
{
  public:
    NemDriver(){};
    virtual ~NemDriver(){}
    static NemDriver* readJSON(std::string fname);
};


class TransferDriver : public NemDriver
{
  public:
    TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                   std::string ofname);

    TransferDriver(std::string srcmsh, std::string trgmsh, std::string method,
                   std::vector<std::string> arrayNames, std::string ofname);

    static TransferDriver* readJSON(json inputjson);
    static TransferDriver* readJSON(std::string ifname);

    ~TransferDriver();

  private:
    meshUser* source;
    meshUser* target; 
};

class RefineDriver : public NemDriver
{
  public:
    RefineDriver(std::string _mesh, std::string method, std::string arrayName,
                 double dev_mult, bool maxIsmin, double edgescale, std::string ofname);
    RefineDriver(std::string _mesh, std::string method, 
                 double edgescale, std::string ofname);
  
  static RefineDriver* readJSON(json inputjson);
  static RefineDriver* readJSON(std::string ifname);

  ~RefineDriver();

  private:
    meshUser* mesh;

};

class MeshGenDriver : public NemDriver
{
  public:
    
    MeshGenDriver(std::string ifname, std::string meshEngine, std::string ofname);
    MeshGenDriver(std::string ifname, std::string meshEngine, 
                  meshingParams* params, std::string ofname);
    ~MeshGenDriver();

    static MeshGenDriver* readJSON(json inputjson);

  private:
    meshUser* mesh; 
    meshingParams* params;
};

#endif
