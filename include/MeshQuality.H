#ifndef MESHQUALITY_H
#define MESHQUALITY_H

#include "nemosys_export.h"
#include "meshBase.H"

#include <vtkMeshQuality.h>
#include <vtkDoubleArray.h>

#ifdef HAVE_CFMSH
  #include "cfmeshQualityParams.H"
  #include "argList.H"
  #include "polyMeshGenModifier.H"
  #include "meshOptimizer.H"
#endif

#include <memory>

class NEMOSYS_EXPORT MeshQuality
{
  public:
    MeshQuality() = default;
    explicit MeshQuality(const meshBase *_mesh);
#ifdef HAVE_CFMSH
    MeshQuality(cfmshQualityParams *params) { _cfmQPrms = params; }
#endif
    ~MeshQuality();

    // disable copy and copy-assignment constructors
    MeshQuality(const MeshQuality &that) = delete;
    MeshQuality &operator=(const MeshQuality &that) = delete;

  public:
    void checkMesh(std::ostream &outputStream);
    void checkMesh();
    void checkMesh(const std::string &fname);
    vtkSmartPointer<vtkDoubleArray> getStats(int n);

  public:
    void cfmOptimize();

  private:
    const meshBase *mesh;
    vtkSmartPointer<vtkMeshQuality> qualityFilter;
#ifdef HAVE_CFMSH
    cfmshQualityParams *_cfmQPrms;
#endif
};

#endif
