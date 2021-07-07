#ifndef NEMOSYS_INPGEOMESH_H_
#define NEMOSYS_INPGEOMESH_H_

#include "Mesh/geoMeshBase.H"

#include <map>
#include <set>
#include <utility>

namespace NEM {
namespace MSH {

/**
 * @brief Class representing meshes in CalculiX input deck (similar to ABAQUS)
 * @details Currently supports the keywords ELEMENT ("C3D4", "C3D8", and "C3D6"
 * types, ELSET parameter), ELSET, NSET, and SURFACE (TYPE=ELEMENT only). By
 * default, physical groups with dimension equal to mesh dimension are turned
 * into ELSETs, and NSETs and physical groups with dimension one less than mesh
 * dimension are turned into SURFACEs and NSETs.
 */
class NEMOSYS_EXPORT inpGeoMesh : public geoMeshBase {
 public:
  static inpGeoMesh *New();
  static inpGeoMesh *Read(const std::string &fileName);
  vtkTypeMacro(inpGeoMesh, geoMeshBase)

  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;

 protected:
  inpGeoMesh();
  explicit inpGeoMesh(const std::string &fileName);

  void reconstructGeo() override;

 private:
  struct NEMOSYS_NO_EXPORT InpSets {
    /**
     * @brief NSET keyword; indexing matches GeoMesh::mesh
     */
    std::map<std::string, std::set<vtkIdType>> nodeSets;
    /**
     * @brief ELSET keyword; indexing matches GeoMesh::mesh
     */
    std::map<std::string, std::set<vtkIdType>> elSets;
    /**
     * @brief SURFACE keyword; indexing matches GeoMesh::sideSet
     */
    std::map<std::string, std::set<vtkIdType>> surfaces;
  };
  /**
   * @brief Ctor meant for ease of delegating to geoMeshBase ctor
   * @param mesh from parsing a .inp mesh file
   */
  explicit inpGeoMesh(std::pair<GeoMesh, InpSets> mesh);
  void resetNative() override;
  static std::pair<GeoMesh, InpSets> inp2GM(const std::string &fileName);

  /**
   * @brief Holds data specific to inp format
   */
  InpSets inpSets_;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_INPGEOMESH_H_
