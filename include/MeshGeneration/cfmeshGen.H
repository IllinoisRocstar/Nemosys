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
#ifndef NEMOSYS_CFMESHGEN_H_
#define NEMOSYS_CFMESHGEN_H_

#ifdef HAVE_CFMSH

#include "MeshGeneration/meshGen.H"

// OpenFOAM headers
#include <argList.H>
#include <fvOptions.H>

class cfmeshParams;

class cfmeshGen : public meshGen {
 public:
  cfmeshGen();
  explicit cfmeshGen(cfmeshParams *params);
  ~cfmeshGen() override;

  // cfmesh mesh creation
 public:
  // create a 3D polygonal discretization from STL file
  int createMeshFromSTL(const char *fname) override;

  // internal management
 private:
  // pre-processor
  void initialize();

  // create dictionaries needed for openFoam
  void createMshDict(const bool &write);

  // surface edge extraction
  int surfaceFeatureEdgeDetect();

  // mesh quality improvement
  int improveMeshQuality();

  // class data and management flags
 private:
  bool defaults;
  std::string caseName{"_cfmesh"};  // internal information
  cfmeshParams *params_;

  // openfoam data structure
 private:
  std::unique_ptr<Foam::dictionary> controlDict_;
  std::unique_ptr<Foam::dictionary> fvSchemes_;
  std::unique_ptr<Foam::dictionary> fvSolution_;
  std::unique_ptr<Foam::Time> runTime_;
  std::unique_ptr<Foam::fvMesh> fmesh_;
  std::unique_ptr<Foam::dictionary> meshDict_;
};

#endif

#endif  // NEMOSYS_CFMESHGEN_H_
