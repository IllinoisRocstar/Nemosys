#ifndef PACKMESHDRIVER_H
#define PACKMESHDRIVER_H

#include "Drivers/NemDriver.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

namespace NEM {
namespace DRV {

/**
  @brief This class drives the flow of pack meshing by using parameters provided
  by users. It supports two different meshing pipelines. First being hexahedral
  conformal meshing and second is tetrahedral periodic as well as non-priodic
  meshing.

  Tetrahedral periodic meshing pipeline operates within rocPack class and uses
  Netgen mesh optimization for final mesh quality improvement if needed.
  ## Pack Meshing Workflow for Tetrahedral Mesh
  > * Geometry generation using RocPack output file.
  > * Enforcing periodicity in geometry for periodic mesh
  > * Mapping periodic surfaces and points with their counterparts.
  > * Meshing the geometry tetrahedrally while enforcing periodic boundary
  boundary conditions.
  > * Exporting final mesh with or without optimization.

**/
class NEMOSYS_EXPORT PackMeshDriver : public NemDriver {
  JSONCONS_TYPE_TRAITS_FRIEND
 protected:
  PackMeshDriver() = default;

 private:
  static constexpr const char *programType = "Pack Mesh Generation";
  jsoncons::string_view getProgramType() const override;
};

}  // namespace DRV
}  // namespace NEM

#endif
