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

class netgenParams;

class netgenGen : public meshGen
{

  public:
    netgenGen();

    netgenGen(netgenParams* params);
    
    ~netgenGen();

    void set_mp(netgenParams* params);
      

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
