#ifndef ICOSIDODECAHEDRONSHAPE_H
#define ICOSIDODECAHEDRONSHAPE_H

#include "nemosys_export.h"
#include <rocPackShape.H>
#include <string>
#include <vector>

namespace NEM {
    
namespace GEO {

/** @brief This class wraps around abstract class rocPackShape. This particular
           shape is called "icosidodecahedron"
**/
class NEMOSYS_EXPORT icosidodecahedronShape : public rocPackShape
{
  // Constructor, Destructor
  public:
  /** @brief icosidodecahedronShape class standard constructor
  **/
  icosidodecahedronShape();

  /** @brief icosidodecahedronShape class standard destructor
  **/
  ~icosidodecahedronShape();

  // Shape Data Methods
  public:
  /** @brief Gets vertices of pack shapes
  **/
  std::vector<std::vector<double>> getVertices() override;

  /** @brief Gets faces of pack shapes
  **/
  std::vector<std::vector<int>> getFaces() override;
  
};
}
}
#endif