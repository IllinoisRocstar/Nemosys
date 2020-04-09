#ifndef NEMOSYS_OSHGEOMESH_H_
#define NEMOSYS_OSHGEOMESH_H_

#include "geoMeshBase.H"

#include "Omega_h_mesh.hpp"

namespace NEM {
namespace MSH {

class NEMOSYS_EXPORT oshGeoMesh : public geoMeshBase {
 public:
  static oshGeoMesh *Read(const std::string &fileName);
  static oshGeoMesh *Read(const std::string &fileName, Omega_h::Library *lib);

 public:
  oshGeoMesh();

  explicit oshGeoMesh(Omega_h::Mesh *oshMesh);
  oshGeoMesh(Omega_h::Mesh *oshMesh, Omega_h::Library *lib);

  ~oshGeoMesh() override;

 public:
  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

 public:
  static vtkUnstructuredGrid *osh2vtk(Omega_h::Mesh *oshMesh);
  static Omega_h::Mesh *vtk2osh(vtkUnstructuredGrid *vtkMesh,
                                Omega_h::Library *lib);

 private:
  Omega_h::Library _oshLibrary;
  Omega_h::Mesh _oshMesh;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_OSHGEOMESH_H_
