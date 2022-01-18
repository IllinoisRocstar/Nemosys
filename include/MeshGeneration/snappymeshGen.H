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
#ifndef NEMOSYS_SNAPPYMESHGEN_H_
#define NEMOSYS_SNAPPYMESHGEN_H_

#ifdef HAVE_CFMSH 

#include "MeshGeneration/meshGen.H"

#include <argList.H>
#include <fvOptions.H>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

class snappymeshParams;

/**
    @brief snappymeshGen facilitates full-hexahedral/hex-dominent meshing of
          complex geometries with surface, patch, or region refinement
          capabilities in addition with adding layers to geometry. A
          background mesh is needed in order to perform meshing operation.
          Inputs can be given in form of surface files (.stl, .obj, etc.).
          Output mesh is written in OpenFOAM polyMesh format.
          snappymeshGen is based on snappyHexMesh utility of OpenFOAM with
          following features
          - Castallated Mesh
          - Snap
          - Add Layers
**/
class snappymeshGen : public meshGen {
  // --- constructors and destructor
 public:
  /** @brief snappymeshGen standard constructor
   **/
  snappymeshGen();  // Default constructor

  /** @brief snappymeshGen alternate constructor
      @param params snappymeshParams object
  **/
  snappymeshGen(snappymeshParams *params);  // Constructors w/ Parameters

  /** @brief snappymeshGen standard destructor
   **/
  ~snappymeshGen();  // destructor

  // --- access
 public:
  /** @brief Creates mesh from input STL file
   **/
  int createMeshFromSTL(const char *fname) override;

 private:
  /** @brief Initializes OpenFOAM environment
   **/
  void initialize();

  /** @brief Creates snappyHexMeshDict for meshing operation
   **/
  void createSnappyDict(const bool &write);

 private:
  /** @brief If enabled, generated mesh using default parameters
   **/
  bool defaults;  // Default scenario

  /** @brief snappymeshParams object Parameters
   **/
  snappymeshParams *params_;  // SnappyHexMesh Parameters

  // -- OpenFOAM Data Structure
 private:
  std::unique_ptr<Foam::dictionary> controlDict_;
  std::unique_ptr<Foam::dictionary> fvSchemes_;
  std::unique_ptr<Foam::dictionary> fvSolution_;
  std::unique_ptr<Foam::dictionary> snappyMshDict_;
  std::unique_ptr<Foam::Time> runTime_;
  std::unique_ptr<Foam::fvMesh> fmesh_;
};

#endif

#endif  // NEMOSYS_SNAPPYMESHGEN_H_
