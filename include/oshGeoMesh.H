#ifndef NEMOSYS_OSHGEOMESH_H_
#define NEMOSYS_OSHGEOMESH_H_

#include "geoMeshBase.H"

// Use forward declaration to avoid importing Omega_h headers
namespace Omega_h {
class Library;
class Mesh;
}  // namespace Omega_h

namespace NEM {
namespace MSH {

/**
 * @class OmegaHInterface
 * @brief management class for @c Omega_h::Library
 * @details To avoid initializing Kokkos multiple times, use this @c
 * OmegaHInterface::GetLibrary() to get an @c Omega_h::Library object
 */
class NEMOSYS_EXPORT OmegaHInterface {
 public:
  ~OmegaHInterface() = default;

  /**@{*/
  OmegaHInterface(const OmegaHInterface &) = delete;
  OmegaHInterface &operator=(const OmegaHInterface &) = delete;
  OmegaHInterface(OmegaHInterface &&) = delete;
  OmegaHInterface &operator=(OmegaHInterface &&) = delete;
  /**@}*/

 public:
  /**
   * Initialize Omega_h. Will only create an Omega_h::Library on first call to
   * method.
   */
  static std::shared_ptr<Omega_h::Library> GetLibrary();

 private:
  OmegaHInterface();

 private:
  static std::shared_ptr<OmegaHInterface> instance;
  std::shared_ptr<Omega_h::Library> library;
};

/**
 * @class oshGeoMesh
 * A concrete implementation of @c geoMeshBase representing a @c Omega_h::Mesh
 */
class NEMOSYS_EXPORT oshGeoMesh : public geoMeshBase {
 public:
  static oshGeoMesh *New();
  vtkTypeMacro(oshGeoMesh, geoMeshBase)

 public:
  /**
   * Create a oshGeoMesh from a file.
   * @param fileName name of file that Omega_h can read
   * @param lib library pointer; if @c nullptr, @c OmegaHInterface::GetLibrary()
   * is used
   * @return new oshGeoMesh
   */
  static oshGeoMesh *Read(const std::string &fileName,
                          Omega_h::Library *lib = nullptr);

 public:
  /**
   * Construct a oshGeoMesh with an empty mesh
   */
  oshGeoMesh();

  /**
   * Construct a oshGeoMesh from an existing Omega_h::Mesh
   * @param oshMesh existing Omega_h::Mesh; note the mesh is copied
   * @param lib library pointer; if @c nullptr, @c OmegaHInterface::GetLibrary()
   * is used
   */
  explicit oshGeoMesh(Omega_h::Mesh *oshMesh, Omega_h::Library *lib = nullptr);

  ~oshGeoMesh() override;

 public:
  void write(const std::string &fileName) override;
  void report(std::ostream &out) const override;
  void takeGeoMesh(geoMeshBase *otherGeoMesh) override;
  void reconstructGeo() override;

 private:
  static GeoMesh osh2GM(Omega_h::Mesh *oshMesh, const std::string &geo = "",
                        const std::string &link = "");
  static Omega_h::Mesh *GM2osh(const GeoMesh &geoMesh,
                               Omega_h::Library *lib = nullptr);

 public:
  /**
   * Get the @c Omega_h::Mesh
   * @return reference to the @c Omega_h::Mesh
   */
  const Omega_h::Mesh &getOshMesh() const { return *_oshMesh; }
  /**
   * Set the mesh to an existing @c Omega_h::Mesh
   * @param oshMesh an existing @c Omega_h::Mesh; note the mesh is actually
   * copied
   */
  void setOshMesh(Omega_h::Mesh *oshMesh);

 private:
  void resetNative() override;

 private:
  std::unique_ptr<Omega_h::Mesh> _oshMesh;
};

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_OSHGEOMESH_H_
