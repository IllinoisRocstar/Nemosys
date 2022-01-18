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
#ifndef NEMOSYS_BLOCKMESHGEN_H_
#define NEMOSYS_BLOCKMESHGEN_H_

#ifdef HAVE_CFMSH

// Nemosys Headers
#include "MeshGeneration/meshGen.H"
#include "Mesh/geoMeshBase.H"
#include "Mesh/foamGeoMesh.H"

// Foam headers
#include <argList.H>
#include <fvOptions.H>

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

class blockMeshParams;  // Friend Class

/**
  @brief blockMeshGen <-- meshGen <-- meshBase
        This class incorporates mesh generating method of blockMesh utility.
        Currently, we are supporting automatic full-hex mesh generation
        Box, Sphere, Cylinder, and Tapered Cone geometries. Output mesh is
        written in OpenFOAM polyMesh format.
**/
class blockMeshGen : public meshGen {
  // --- constructors and destructor
 public:
  /** @brief blockMeshGen standard constructor
   **/
  blockMeshGen();  // Default constructor

  /** @brief blockMeshGen alternate constructor. Uses user-defined parameters
            to perform requested meshing operation.
      @param params blockMeshParams object
  **/
  blockMeshGen(blockMeshParams *params);  // Constructor with parameters

  /** @brief blockMeshGen standard desctructor
   **/
  ~blockMeshGen();  // Class destructor

  // -- meshGen methods override
 public:
  /** @brief Generates mesh and returns VTK database
      @param fname Input surface file name
      @return meshBase mesh dataset
  **/
  int createMeshFromSTL(const char *fname);

  // -- Internal
 private:
  /**
      @brief Initializes OpenFOAM args and runtime
  **/
  void initialize();

  /**
      @brief Creates blockMeshDict from user arguments
  **/
  void createBlockMshDict(const bool &write);

  // --- Internal Use
 private:
  bool defaults;  // Default parameters boolean

  /** @brief Definition of parameters pointer
   **/
  blockMeshParams *params_;  // Definition of parameters pointer

  // --- OpenFOAM Data Structure
 private:
  std::unique_ptr<Foam::Time> runTime_;
  std::unique_ptr<Foam::fvMesh> fmesh_;
  std::unique_ptr<Foam::dictionary> blockMshDict_;
  std::unique_ptr<Foam::dictionary> controlDict_;
  std::unique_ptr<Foam::dictionary> fvSchemes_;
  std::unique_ptr<Foam::dictionary> fvSolution_;
};

#endif

#endif  // NEMOSYS_BLOCKMESHGEN_H_
