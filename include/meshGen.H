#ifndef MESHGEN_H
#define MESHGEN_H

#include <meshBase.H>

class meshingParams
{
  public:
    meshingParams(){};
    virtual ~meshingParams(){};    
    
};

class NetgenParams : public meshingParams
{
  
  public:
    // initializing params with defaults
    NetgenParams();
    ~NetgenParams(){};
  
  public:
    bool uselocalh;
    double maxh;
    double fineness;  
    double grading;
    double elementsperedge;
    double elementspercurve;
    bool closeedgeenable;
    double closeedgefact;
    bool second_order;
    std::string meshsize_filename;
    bool quad_dominated;
    bool optvolmeshenable;
    int optsteps_2d;
    int optsteps_3d;
    bool invert_tets;
    bool invert_trigs;
    bool check_overlap;
    bool check_overlapping_boundary;
};

class meshGen 
{
  public:
    meshGen(){};
    virtual ~meshGen(){}
    
    // creates generator with default parameters
    static meshGen* Create(std::string fname, std::string meshEngine);
    // creates generater with specified parameters
    static meshGen* Create(std::string fname, std::string meshEngine, meshingParams* params);
    virtual int createMeshFromSTL(char* fname) = 0;

};

class meshNetgen : public meshGen
{

  public:
    meshNetgen();

    meshNetgen(NetgenParams* params);
    
    ~meshNetgen();

    void set_mp(NetgenParams* params);
      

  // netgen mesh creation
  public:
    int createMeshFromSTL(char* fname);
  private:
    nglib::Ng_Meshing_Parameters mp; // params for netgen meshing  
    nglib::Ng_Mesh* mesh; // netgen mesh object
};

#endif
