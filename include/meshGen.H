#ifndef MESHGEN_H
#define MESHGEN_H

#include <meshBase.H>

class meshGen 
{
  public:
    meshGen(){};
    ~meshGen(){}; 
    virtual int createMeshFromSTL(char* fname) = 0;
};

class meshNetgen : public meshGen
{
  // constructor and destructor
  public:
    meshNetgen():mp()
    {
      std::cout << "initializing netgen mesh generator" << std::endl;
      nglib::Ng_Init();
      mesh = nglib::Ng_NewMesh();
    }
    
    ~meshNetgen()
    {
      std::cout << "finalizing netgen mesh generator" << std::endl;
      if(mesh) nglib::Ng_DeleteMesh(mesh);
      nglib::Ng_Exit();
    }      

  // netgen mesh creation
  public:
    int createMeshFromSTL(char* fname);
  private:
    nglib::Ng_Meshing_Parameters mp; // params for netgen meshing  
    nglib::Ng_Mesh* mesh; // netgen mesh object
};

#endif
