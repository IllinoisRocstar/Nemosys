#include <iostream>
#include <rocPackShape.H>
#include <hmxShape.H>
#include <petnShape.H>
#include <icosidodecahedronShape.H>
#include <memory>


namespace NEM {

namespace GEO {

rocPackShape *rocPackShape::getShape(const std::string &shapeName)
{
  if (shapeName == "hmx")
  {
    hmxShape *assignShape = new hmxShape();
    return assignShape;
   }
  else if (shapeName == "petn")
  {
    petnShape *assignShape = new petnShape();
    return assignShape;
  }
  else if (shapeName == "icosidodecahedron")
  {
    icosidodecahedronShape *assignShape = new icosidodecahedronShape();
    return assignShape;
  }
  else
  {
    std::cerr << "The " << shapeName << " shape is not supported yet!"
              << std::endl;
    throw;
  }
}

}

}