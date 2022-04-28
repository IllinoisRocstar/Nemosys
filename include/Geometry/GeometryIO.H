#ifndef NEMOSYS_GEOMETRYIO_H_
#define NEMOSYS_GEOMETRYIO_H_

#include "nemosys_export.h"

#include "Geometry/GeoManager.H"

// class TopoDS_Shape;
// class BRepBuilderAPI_Sewing;

namespace NEM {
namespace GEO {

/**
 * List of supported @c GeometryIOBase classes
 *
 * Note: When adding a new format, extend the switch statements in the static
 * factory methods.
 */
enum class NEMOSYS_EXPORT GeomType { STEP_GEO, IGES_GEO, STL_GEO };

/**
 * @brief Class to manage @c TopoDS_Shapes along with metadata
 * @details Handles results of @c BRepBuilderAPI_MakeShape (in particular,
 * boolean operations) by delegating to the arguments' metadata. Also see @c
 * NEM::GEO::ShapeData
 */
class NEMOSYS_EXPORT GeometryIO {
 public:
  /**
   * Create an empty geometry manager
   * @param dim Shapes of this dimension (clipped to be between 1 and 3,
   * inclusive) are treated as first-class citizens, with lower-dimensional
   * shapes assumed to be subshapes.
   */
  explicit GeometryIO(GeoManager *geoManager, std::string &filename);

  ~GeometryIO();

  TopoDS_Shape sewShapes();

  void write(std::string &filename);

  void writeSTEP(std::string &filename, TopoDS_Shape &sewnShape);

  void writeIGES(std::string &filename, TopoDS_Shape &sewnShape);

  void writeSTL(std::string &filename, TopoDS_Shape &sewnShape);

  // bool isChild(const TopoDS_Shape &shape) const;

 private:
  GeoManager *geoManager_;
};

}  // namespace GEO
}  // namespace NEM

#endif  // NEMOSYS_GEOMETRYIO_H_
