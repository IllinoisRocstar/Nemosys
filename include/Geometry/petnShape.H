#ifndef NEMOSYS_PETNSHAPE_H_
#define NEMOSYS_PETNSHAPE_H_

#include "nemosys_export.h"
#include "Geometry/rocPackShape.H"
#include <string>
#include <vector>

namespace NEM {
    
namespace GEO {

/** @brief This class wraps around abstract class rocPackShape. This particular
           shape is called "PETN", which is a shape of "Pentaerythritol 
           tetranitrate" crystal.
           Ref -> https://en.wikipedia.org/wiki/Pentaerythritol_tetranitrate
**/
class NEMOSYS_EXPORT petnShape : public rocPackShape
{
  // Constructor, Desctructor
  public:
  /** @brief petnShape class standard constructor
  **/
  petnShape();

  /** @brief petnShape class standard destructor
  **/
  ~petnShape();

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

#endif  // NEMOSYS_PETNSHAPE_H_
