#ifndef NEMOSYS_GMSHGEOMESH_H_
#define NEMOSYS_GMSHGEOMESH_H_

#include "geoMeshBase.H"

namespace NEM {
namespace MSH {

class NEMOSYS_EXPORT gmshGeoMesh : public geoMeshBase {
 public:
  static gmshGeoMesh *Read(const std::string &fileName);

 public:
  gmshGeoMesh();

  explicit gmshGeoMesh(const std::string &gmshMesh);

  ~gmshGeoMesh() override;

 public:
  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

 private:
  static GeoMesh gmsh2GM(const std::string &gmshMesh);
  static std::string GM2gmsh(const GeoMesh &geoMesh);

 private:
  std::string _gmshMesh;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_GMSHGEOMESH_H_
