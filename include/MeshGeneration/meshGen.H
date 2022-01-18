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
#ifndef NEMOSYS_MESHGEN_H_
#define NEMOSYS_MESHGEN_H_

#include "nemosys_export.h"
#include "Mesh/geoMeshFactory.H"
#include "Mesh/geoMeshBase.H"

#include <vtkSmartPointer.h>
#include <vtkDataSet.h>

#include <string>
#include <iostream>

class NEMOSYS_EXPORT meshGen {
 public:
  meshGen() : dataSet(nullptr) {
    std::cout << "Base mesh generator created" << std::endl;
  }

  virtual ~meshGen() {
    std::cout << "Base mesh generator destroyed" << std::endl;
  }

  virtual int createMeshFromSTL(const char *fname) = 0;
  vtkSmartPointer<vtkDataSet> getDataSet() const { return dataSet; }
  std::unique_ptr<NEM::MSH::geoMeshBase> gmData;

 protected:
  vtkSmartPointer<vtkDataSet> dataSet;
};

#endif  // NEMOSYS_MESHGEN_H_
