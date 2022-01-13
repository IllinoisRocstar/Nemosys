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
