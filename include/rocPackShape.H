#ifndef ROCPACKSHAPE_H
#define ROCPACKSHAPE_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "nemosys_export.h"

namespace NEM {

namespace GEO {

/** @brief This is an abstract class for all shapes supported by RocPack
         utility. This class contains abstract methods used for defining
         various shapes. All methods are defined in various shape classes
         under NEM::GEO namespace.

         # The shapes currently supported are,
         > * HMX (Octogen crystal)
         > * PETN (Pentaerythritol tetranitrate crystal)
         > * ICOSIDODECAHEDRON
**/
class NEMOSYS_EXPORT rocPackShape {
  // Constructor & Deconstructor
 public:
  /** @brief Default constructor
   **/
  rocPackShape() {}

  /** @brief Standard Deconstructor
   **/
  virtual ~rocPackShape() {
    std::cout << "rocPackShape class destroyed!" << std::endl;
  }

  // Access
  // public:
  /** @brief Creates shape object for requsted shape
      @param shapeName Name of shape requested
      @retun Shape pointer
  **/
  static std::shared_ptr<rocPackShape> getShape(const std::string &shapeName);

  // Access
 public:
  /** @brief Gets vertices of pack shapes
   **/
  virtual std::vector<std::vector<double>> getVertices() = 0;

  /** @brief Gets faces of pack shapes
   **/
  virtual std::vector<std::vector<int>> getFaces() = 0;
};

}  // namespace GEO

}  // namespace NEM

#endif
