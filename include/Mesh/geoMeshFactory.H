#ifndef NEMOSYS_GEOMESHFACTORY_H_
#define NEMOSYS_GEOMESHFACTORY_H_

#include "Mesh/geoMeshBase.H"

namespace NEM {
namespace MSH {

/**
 * List of supported @c geoMeshBase classes
 *
 * Note: When adding a new format, extend the switch statements in the static
 * factory methods.
 */
enum class NEMOSYS_EXPORT MeshType {
  VTK_GEO_MESH,  ///< Based on vtkUnstructuredGrid from VTK library.
  GMSH_GEO_MESH, ///< Mesh and geometry based on Gmsh's public API.
  OSH_GEO_MESH,  ///< Based on Omega_h::Mesh from Omega_h library.
  EXO_GEO_MESH,  ///< Interface to EXODUS II library.
  FOAM_GEO_MESH, ///< Based on fvMesh from OpenFOAM library.
  INP_GEO_MESH,  ///< Abaqus/CalculiX mesh input format
  SMESH_GEO_MESH,///< Based on SMESH_Mesh from Salome SMESH (see contrib/)
};

/**
 * Read a mesh from file
 * @param fileName name of file
 * @param type type of mesh
 * @return pointer to base mesh class
 */
NEMOSYS_EXPORT geoMeshBase *Read(const std::string &fileName, MeshType type);
/**
 * Read a mesh from file. Determines mesh type from file extension.
 * @param fileName name of file
 * @return pointer to base mesh class
 */
NEMOSYS_EXPORT geoMeshBase *Read(const std::string &fileName);

/**
 * Create a new mesh object
 * @param type type of mesh
 * @return pointer to base mesh class
 */
NEMOSYS_EXPORT geoMeshBase *New(MeshType type);
/**
 * Create a new mesh object. Determines mesh type from file extension.
 * @param fileName name of file
 * @return pointer to base mesh class
 */
NEMOSYS_EXPORT geoMeshBase *New(const std::string &fileName);

}  // namespace MSH
}  // namespace NEM

#endif  // NEMOSYS_GEOMESHFACTORY_H_
