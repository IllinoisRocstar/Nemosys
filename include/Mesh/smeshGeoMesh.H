#ifndef NEMOSYS_SMESHGEOMESH_H_
#define NEMOSYS_SMESHGEOMESH_H_

#include "nemosys_export.h"
#include "Mesh/geoMeshBase.H"

#include <memory>

class SMESH_Gen;
class SMESH_Mesh;

namespace NEM {
namespace MSH {

class NEMOSYS_EXPORT smeshGeoMesh : public geoMeshBase {
 public:
  smeshGeoMesh();
  ~smeshGeoMesh() override;
  static smeshGeoMesh *New();
  vtkTypeMacro(smeshGeoMesh, geoMeshBase)

  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

  void setSMeshMesh(std::unique_ptr<SMESH_Mesh> &&mesh,
                    std::shared_ptr<SMESH_Gen> gen);
  const SMESH_Mesh &getSMESHMesh() const;

 private:
  void resetNative() override;

  static GeoMesh SmeshToGM(SMESH_Mesh &mesh);
  static void GMToSMESH(const GeoMesh& geoMesh, SMESH_Mesh &outMesh);

  std::shared_ptr<SMESH_Gen> gen_;
  std::unique_ptr<SMESH_Mesh> mesh_;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_SMESHGEOMESH_H_
