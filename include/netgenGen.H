#ifndef NETGENGEN_H
#define NETGENGEN_H
#include <meshGen.H>
namespace nglib
{
  #include<nglib.h>
}

// TODO: add in-mem conversion to vtkdataSet and get method.
//       in that way, we can use the 'copy' constructor from meshbase 
//       and don't need to use the file-based export method

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
    bool refine_with_geom;
    bool refine_without_geom;
};

class netgenGen : public meshGen
{

  public:
    netgenGen();

    netgenGen(NetgenParams* params);
    
    ~netgenGen();

    void set_mp(NetgenParams* params);
      

  // netgen mesh creation
  public:
    int createMeshFromSTL(const char* fname);
  private:
    nglib::Ng_Meshing_Parameters mp; // params for netgen meshing  
    nglib::Ng_Mesh* mesh; // netgen mesh object
    bool refine_with_geom; // if refinement enabled, adapt to geom
    bool refine_without_geom; // if refinement enabled, just do uniform
};

#endif
