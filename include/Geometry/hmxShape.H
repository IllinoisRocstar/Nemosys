#ifndef NEMOSYS_HMXSHAPE_H_
#define NEMOSYS_HMXSHAPE_H_

#include "nemosys_export.h"
#include "Geometry/rocPackShape.H"
#include <string>
#include <vector>

namespace NEM {
    
namespace GEO {

/** @brief This class wraps around abstract class rocPackShape. This particular
           shape is called "HMX", which is a shape of "Octogen" crystal.
           Ref -> https://en.wikipedia.org/wiki/HMX
**/
class NEMOSYS_EXPORT hmxShape : public rocPackShape
{
  // Constructor, Destructor
  public:
  /** @brief hmxShape class standard constructor
  **/
  hmxShape();

  /** @brief hmxShape class standard destructor
  **/
  ~hmxShape();

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
#endif  // NEMOSYS_HMXSHAPE_H_
