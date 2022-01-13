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
#ifndef NEMOSYS_NETGENGEN_H_
#define NEMOSYS_NETGENGEN_H_

#include "nemosys_export.h"
#include "MeshGeneration/meshGen.H"

namespace nglib
{
  #include <nglib.h>
}

// TODO: add in-mem conversion to vtkDataSet and get method.
//       in that way, we can use the 'copy' constructor from meshBase
//       and don't need to use the file-based export method

class netgenParams;

class NEMOSYS_EXPORT netgenGen : public meshGen
{
  public:
    netgenGen();
    explicit netgenGen(const netgenParams *params);

    ~netgenGen() override;

    void set_mp(const netgenParams *params);

  // netgen mesh creation
  public:
    int createMeshFromSTL(const char *fname) override;

  private:
    nglib::Ng_Meshing_Parameters mp; // params for netgen meshing  
    nglib::Ng_Mesh* mesh; // netgen mesh object
    bool refine_with_geom; // if refinement enabled, adapt to geom
    bool refine_without_geom; // if refinement enabled, just do uniform
};

#endif  // NEMOSYS_NETGENGEN_H_
