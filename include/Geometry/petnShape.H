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
