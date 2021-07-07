#ifndef NEMOSYS_GMSHGEOMESH_H_
#define NEMOSYS_GMSHGEOMESH_H_

#include "Mesh/geoMeshBase.H"

namespace NEM {
namespace MSH {

/**
 * @class gmshGeoMesh
 * A concrete implementation of @c geoMeshBase representing a Gmsh mesh.
 */
class NEMOSYS_EXPORT gmshGeoMesh : public geoMeshBase {
 public:
  static gmshGeoMesh *New();
  vtkTypeMacro(gmshGeoMesh, geoMeshBase)

 public:
  /**
   * Create a gmshGeoMesh from a file
   * @param fileName file name that gmsh can read
   * @return new gmshGeoMesh object
   */
  static gmshGeoMesh *Read(const std::string &fileName);
  /**
   * Convert gmsh cell type to VTK cell type
   * @param gmshType gmsh cell type
   * @return VTK cell type
   */
  static VTKCellType getVTKTypeFromGmshType(const std::string &gmshType);

 public:
  /**
   * Construct a gmshGeoMesh with an empty mesh.
   */
  gmshGeoMesh();

  /**
   * Construct a gmshGeoMesh from an existing gmsh model
   * @param gmshMesh name of gmsh model created from @c gmsh::open or @c
   * gmsh::add
   */
  explicit gmshGeoMesh(const std::string &gmshMesh);

  ~gmshGeoMesh() override;

 public:
  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

 private:
  /**
   * Convert @p gmshMesh (name of a gmsh model) to GeoMesh
   */
  static GeoMesh gmsh2GM(const std::string &gmshMesh);
  /**
   * Add the mesh from @p geoMesh to the gmsh model in @p geoMesh
   */
  static std::string GM2gmsh(const GeoMesh &geoMesh);

  void resetNative() override;

 private:
  std::string _gmshMesh;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_GMSHGEOMESH_H_
