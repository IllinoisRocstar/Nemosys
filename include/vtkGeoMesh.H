#ifndef NEMOSYS_VTKGEOMESH_H_
#define NEMOSYS_VTKGEOMESH_H_

#include "geoMeshBase.H"

namespace NEM {
namespace MSH {

class NEMOSYS_EXPORT vtkGeoMesh : public geoMeshBase {
 public:
  static vtkGeoMesh *Read(const std::string &fileName,
                          const std::string &phyGrpArrayName = std::string());

 public:
  vtkGeoMesh();

  explicit vtkGeoMesh(vtkUnstructuredGrid *inUnstructuredGrid,
                      const std::string &phyGrpArrayName = std::string());

  vtkGeoMesh(const vtkGeoMesh &) = delete;
  vtkGeoMesh &operator=(const vtkGeoMesh &) = delete;

  vtkGeoMesh(vtkGeoMesh &&) = default;
  vtkGeoMesh &operator=(vtkGeoMesh &&) = default;

  ~vtkGeoMesh() override;

 public:
  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

 private:
  static GeoMesh vtk2GM(vtkUnstructuredGrid *vtkMesh,
                        const std::string &phyGrpArrayName = std::string());
  static vtkUnstructuredGrid *GM2vtk(const GeoMesh &geoMesh);

 private:
  vtkSmartPointer<vtkUnstructuredGrid> _vtkMesh;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_VTKGEOMESH_H_
