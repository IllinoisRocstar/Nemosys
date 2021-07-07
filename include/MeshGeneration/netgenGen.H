#ifndef NEMOSYS_NETGENGEN_H_
#define NEMOSYS_NETGENGEN_H_

#include "nemosys_export.h"
#include "MeshGeneration/meshGen.H"

#ifdef HAVE_NGEN
namespace nglib
{
  #include <nglib.h>
}
#endif

// TODO: add in-mem conversion to vtkDataSet and get method.
//       in that way, we can use the 'copy' constructor from meshBase
//       and don't need to use the file-based export method

class netgenParams;

class NEMOSYS_EXPORT netgenGen : public meshGen
{
  public:
    netgenGen();
    explicit netgenGen(const netgenParams *params);

    ~netgenGen() override;

    void set_mp(const netgenParams *params);

  // netgen mesh creation
  public:
    int createMeshFromSTL(const char *fname) override;

  private:
#ifdef HAVE_NGEN
    nglib::Ng_Meshing_Parameters mp; // params for netgen meshing  
    nglib::Ng_Mesh* mesh; // netgen mesh object
#endif
    bool refine_with_geom; // if refinement enabled, adapt to geom
    bool refine_without_geom; // if refinement enabled, just do uniform
};

#endif  // NEMOSYS_NETGENGEN_H_
