#ifndef NEMOSYS_REFINE_H_
#define NEMOSYS_REFINE_H_

#include <string>

#include <MAdLib.h>

#include "nemosys_export.h"
#include "Mesh/meshBase.H"

namespace NEM {
namespace ADP {

class NEMOSYS_EXPORT Refine {
  // constructors and destructors
 public:
  Refine(meshBase *_mesh, const std::string &method, int arrayID,
         double dev_mult, bool maxIsmin, double edge_scale,
         const std::string &_ofname, double sizeFactor = 1.0, int order = 1);
  ~Refine();

  // methods
 public:
  void run(bool transferData, bool bndryConstraint = false);

 private:
  meshBase *mesh;      // mesh to be refined
  MAd::pMesh MadMesh;  // MAdLib mesh object generated from converted meshBase
  MAd::BackgroundSF *bSF;     // background sizeField
  MAd::PWLSField *pwlSF;      // piecewise linear size field
  MAd::MeshAdapter *adapter;  // adapter
  std::string ofname;
  MAd::pGEntity bnd;
  bool bndrConst;

  // helpers
 private:
  void initUniform(double edge_scale);
  void initAdaptive(int arrayID, const std::string &method);

  void classifyBoundaries();
  void unClassifyBoundaries();
};

}  // namespace ADP
}  // namespace NEM

#endif  // NEMOSYS_REFINE_H_
