#ifndef PETNSHAPE_H
#define PETNSHAPE_H

#include "nemosys_export.h"
#include <rocPackShape.H>
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
#endif