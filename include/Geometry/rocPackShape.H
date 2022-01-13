/*******************************************************************************
* Promesh                                                                      *
* Copyright (C) 2022, IllinoisRocstar LLC. All rights reserved.                *
*                                                                              *
* Promesh is the property of IllinoisRocstar LLC.                              *
*                                                                              *
* IllinoisRocstar LLC                                                          *
* Champaign, IL                                                                *
* www.illinoisrocstar.com                                                      *
* promesh@illinoisrocstar.com                                                  *
*******************************************************************************/
/*******************************************************************************
* This file is part of Promesh                                                 *
*                                                                              *
* This version of Promesh is free software: you can redistribute it and/or     *
* modify it under the terms of the GNU Lesser General Public License as        *
* published by the Free Software Foundation, either version 3 of the License,  *
* or (at your option) any later version.                                       *
*                                                                              *
* Promesh is distributed in the hope that it will be useful, but WITHOUT ANY   *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    *
* FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more *
* details.                                                                     *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this program. If not, see <https://www.gnu.org/licenses/>.        *
*                                                                              *
*******************************************************************************/
#ifndef NEMOSYS_ROCPACKSHAPE_H_
#define NEMOSYS_ROCPACKSHAPE_H_

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

#endif  // NEMOSYS_ROCPACKSHAPE_H_
