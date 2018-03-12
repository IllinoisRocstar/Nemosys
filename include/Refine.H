#ifndef REFINE_H
#define REFINE_H

#include <meshBase.H>

class Refine
{
  // constructors and destructors
  public:
    Refine(){}
    Refine(meshBase* _mesh, const std::string& method, 
           int arrayID, double dev_mult, bool maxIsmin, 
           double edge_scale, std::string _ofname);
    ~Refine();

  // methods
  public:
    void run(bool transferData);

  private:
    meshBase* mesh; // mesh to be refined
    MAd::pMesh MadMesh; // MAdLib mesh object generated from converted meshBase 
    MAd::BackgroundSF* bSF; // background SizeField 
    MAd::PWLSField* pwlSF; // piecewise linear size field
    MAd::MeshAdapter* adapter; // adapter
    std::string ofname;
 
  // helpers
  private:
    void initUniform(double edge_scale);
    void initAdaptive(int arrayID, const std::string& method);
  
    void classifyBoundaries();
    void unClassifyBoundaries();
};

#endif
