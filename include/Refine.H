#ifndef REFINE_H
#define REFINE_H

#include "nemosys_export.h"
#include "meshBase.H"

// MAdLib
#include <MAdLib.h>

#include <string>

class NEMOSYS_EXPORT Refine
{
  // constructors and destructors
  public:
    Refine() = default;

    Refine(meshBase *_mesh, const std::string &method,
           int arrayID, double dev_mult, bool maxIsmin, double edge_scale,
           const std::string &_ofname,
           double sizeFactor = 1.);
    ~Refine();

  // methods
  public:
    void run(bool transferData);

  private:
    meshBase *mesh; // mesh to be refined
    MAd::pMesh MadMesh; // MAdLib mesh object generated from converted meshBase 
    MAd::BackgroundSF *bSF; // background sizeField
    MAd::PWLSField *pwlSF; // piecewise linear size field
    MAd::MeshAdapter *adapter; // adapter
    std::string ofname;

  // helpers
  private:
    void initUniform(double edge_scale);
    void initAdaptive(int arrayID, const std::string &method);

    void classifyBoundaries();
    void unClassifyBoundaries();
};

#endif
