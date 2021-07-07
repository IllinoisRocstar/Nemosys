#include "Geometry/hmxShape.H"
#include "Geometry/icosidodecahedronShape.H"
#include "Geometry/petnShape.H"
#include "Geometry/rocPackShape.H"

namespace NEM {

namespace GEO {

std::shared_ptr<rocPackShape> rocPackShape::getShape(
    const std::string &shapeName) {
  if (shapeName == "hmx") {
    std::shared_ptr<hmxShape> assignShape(new hmxShape());
    return assignShape;
  } else if (shapeName == "petn") {
    std::shared_ptr<petnShape> assignShape(new petnShape());
    return assignShape;
  } else if (shapeName == "icosidodecahedron") {
    std::shared_ptr<icosidodecahedronShape> assignShape(
        new icosidodecahedronShape());
    return assignShape;
  } else {
    std::cerr << "The " << shapeName << " shape is not supported yet!"
              << std::endl;
    throw;
  }
}

}  // namespace GEO

}  // namespace NEM
