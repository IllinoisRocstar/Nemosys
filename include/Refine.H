#ifndef REFINE_H
#define REFINE_H

#include <meshBase.H>

class Refine
{
  // constructors and destructors
  public:
    Refine(){}
    Refine(meshBase* _mesh, std::string method, 
           int arrayID, double dev_mult, bool maxIsmin);
    ~Refine();

  // methods
  public:
    void run();   


  private:
    meshBase* mesh;
    MAd::pMesh MadMesh;
    MAd::BackgroundSF* bSF;
    MAd::MeshAdapter* adapter;
  
  // helpers
  private:
    void classifyBoundaries();
    void unClassifyBoundaries();
};

#endif
