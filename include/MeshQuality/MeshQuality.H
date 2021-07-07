#ifndef NEMOSYS_MESHQUALITY_H_
#define NEMOSYS_MESHQUALITY_H_

#include "nemosys_export.h"

#include <ostream>
#include <string>

#include <vtkDoubleArray.h>
#include <vtkMeshQuality.h>
#include <vtkSmartPointer.h>
#ifdef HAVE_CFMSH
#include "MeshQuality/cfmeshQualityParams.H"
#endif

#include "Mesh/meshBase.H"

class NEMOSYS_EXPORT MeshQuality {
 public:
  MeshQuality() = default;
  explicit MeshQuality(const meshBase *_mesh);
#ifdef HAVE_CFMSH
  explicit MeshQuality(const cfmeshQualityParams *params) { _cfmQPrms = params; }
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
  const cfmeshQualityParams *_cfmQPrms;
#endif
};

#endif  // NEMOSYS_MESHQUALITY_H_
