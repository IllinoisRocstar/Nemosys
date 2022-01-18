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
#ifndef NEMOSYS_GMSHGEN_H_
#define NEMOSYS_GMSHGEN_H_

#include "MeshGeneration/meshGen.H"

namespace NEM {

namespace GEN {

class gmshParams;

class gmshGen : public meshGen {
 public:
  /** @brief gmshGen default constructor
   **/
  gmshGen();

  /** @brief gmshGen alternate constructor with parameters
@param params gmshParams object
**/
  gmshGen(gmshParams *params);

  /** @brief gmshGen standard destructor
   **/
  ~gmshGen();

 public:
  /** @brief Creates mesh from input STEP file
      @param fname The input STEP file
  **/
  int createMeshFromSTL(const char *fname) override;

  int createMeshFromSTEP(const char *fname);

 public:

  /** @brief Get filename of generated .msh file
   **/
  const std::string &getMshFilename() const { return mshFname; }

 private:
  /** @brief Sets the global geometry and meshing options
   **/
  void globalOptions();

  /** @brief Gets geometry entitiy names of STEP geometry
   **/
  void getGeomNames();

  /** @brief Gets the surface colors of STEP geometry
   **/
  void getSurfaceColors();

  /** @brief Applies mesh size fields
   **/
  void meshSizeFields();

  /** @brief Applies physical names based on color
   **/
  void applyColorNames();

  /** @brief Applies transfinite settings to prescribed hexahedral volumes
   **/
  void applyTransfiniteVolumes();

  /** @brief gmshParams object Parameters
   **/
  gmshParams *meshParams;  // Gmsh Parameters

  /** @brief Name of .msh file output by gmsh
   */
  std::string mshFname;
};

}  // namespace GEN

}  // namespace NEM

#endif  // NEMOSYS_GMSHGEN_H_
