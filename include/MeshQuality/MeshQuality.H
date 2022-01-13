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
#ifndef NEMOSYS_MESHQUALITY_H_
#define NEMOSYS_MESHQUALITY_H_

#include "nemosys_export.h"

#include <ostream>
#include <string>

#include <vtkDoubleArray.h>
#include <vtkMeshQuality.h>
#include <vtkSmartPointer.h>
#ifdef HAVE_CFMSH
#include "MeshQuality/cfmeshQualityParams.H"
#endif

#include "Mesh/meshBase.H"

class NEMOSYS_EXPORT MeshQuality {
 public:
  MeshQuality() = default;
  explicit MeshQuality(const meshBase *_mesh);
#ifdef HAVE_CFMSH
  explicit MeshQuality(const cfmeshQualityParams *params) { _cfmQPrms = params; }
#endif
  ~MeshQuality();

  // disable copy and copy-assignment constructors
  MeshQuality(const MeshQuality &that) = delete;
  MeshQuality &operator=(const MeshQuality &that) = delete;

 public:
  void checkMesh(std::ostream &outputStream);
  void checkMesh();
  void checkMesh(const std::string &fname);
  vtkSmartPointer<vtkDoubleArray> getStats(int n);

 public:
  void cfmOptimize();

 private:
  const meshBase *mesh;
  vtkSmartPointer<vtkMeshQuality> qualityFilter;
#ifdef HAVE_CFMSH
  const cfmeshQualityParams *_cfmQPrms;
#endif
};

#endif  // NEMOSYS_MESHQUALITY_H_
