#ifndef NEMOSYS_MESHMANIPULATIONFOAM_H_
#define NEMOSYS_MESHMANIPULATIONFOAM_H_

#include "nemosys_export.h"

// VTK
#include <vtkDoubleArray.h>
#include <vtkMeshQuality.h>

// OpenFOAM
#include <argList.H>
#include <fileName.H>
#include <fvMesh.H>

// Nemosys
#include "MeshManipulationFoam/MeshManipulationFoamParams.H"
#include "Mesh/geoMeshFactory.H"
#include "Mesh/geoMeshBase.H"
#include "Mesh/foamGeoMesh.H"
#include "Mesh/vtkGeoMesh.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

/** @brief MeshManipulation class contains several OpenFOAM utilities to perform
    various mesh manipulation operations. These utilities reads Foam mesh
    and outputs Foam mesh. Some of them are surface mesh manipulation
    utilities.
    Manipulation utilities include
    - surfaceLambdaMuSmooth
    - splitMeshRegions
    - mergeMeshes
    - createPatch
    - foamToSurface
    - surfaceSplitByTopology.

    Foam mesh reading and writing capabilities are built inside each
manipulation utility.
**/

class NEMOSYS_EXPORT MeshManipulationFoam {
  // --- constructors and destructors
 public:
  /**
      @brief Standard MeshManipulationFoam Constructor
  **/
  MeshManipulationFoam(){};  // Standard Constructor

  /**
      @brief Custom MeshManipulationFoam Constructor
      @param params MeshManipulationFoamParams Object
  **/
  MeshManipulationFoam(MeshManipulationFoamParams *params) {
    _mshMnipPrms = params;
  }  // Alternate Constructor

  /**
      @brief Standard MeshManipulationFoam Destructor
  **/
  ~MeshManipulationFoam();  // Class Destructor

  // --- meshBase factory methods
 public:
  /** @brief surfLambdaMuSmooth performs laplacian smoothing over surface
            and writes it into output surface file. This utility supports .ofs,
            .obj, .inp, .stl, .tri, .off, .stlb, .nas, .bdf, .gts, .vtk, .ac
            extensions for input and output surface files.
  **/
  void surfLambdaMuSmooth();  // SurfaceLambdaMuSmooth utility

  /** @brief splitMshRegions walks through mesh using cell-face-cell walk and
            identifies different cellZones as different regions. Writes these
            regions in separate Foam meshes.
      @param
      @return domain number skipped during splitting, Total disconnected regions
  **/
  std::pair<std::vector<int>, std::vector<std::string>> splitMshRegions();

  /** @brief mergeMesh method merges two fvMesh datasets into one fvMesh
             including the patches and domain information
      @param dirStat domain number to skip (output from splitMshRegions)
      @param nDomains Number of meshes to merge
  **/
  void mergeMesh(int dirStat, int nDomains);  // MergeMeshes utility

  /** @brief createPtch utility creates user-defined patches in Foam mesh. User
            input is usually createPatchDict.
      @param dirStat domain number to skip (output from splitMshRegions)
  **/
  void createPtch(int dirStat);  // CreatePatch utility

  /** @brief foamToSurf utility reads Foam mesh, extracts its surface and
            writes into an STL file.
  **/
  void foamToSurf();  // FoamToSurface utility

  /** @brief surfSpltByTopology utility takes a surface file as an input and
            automatically divides disconnected regions into separate patches.
            writes a surface file with separated patches as well as multiple
            surface files for separate patches.
  **/
  int surfSpltByTopology();  // surfaceSplitByTopology utility

  /** @brief addCohesiveElements method supports addition of zero-sized
            connectivity elements at conformal shared interface between two
            meshes. This method can take two .vtk/.vtu files or Foam meshes as
  input and outputs a single .vtk/.vtu file with original meshes plus cohesive
  elements.
      @param tol Tolerance for conformal duplicate nodes between two meshes.
      @param outName Output mesh name with extension (i.e "output.vtu")
  **/
  // Adds connectivity info at duplicate nodes.
  void addCohesiveElements(double tol, const std::string &outName);

  /** @brief This method adds connectivity elements with finite thickness
        at shared interface between two conformal meshes
      @param tol Tolerance for conformal duplicate nodes between two meshes.
      @param outName Output mesh name with extension (i.e "output.vtu")
      @param thickness Desired thickness for new cohesive elements.
  **/
  void addArtificialThicknessElements(double &tol, const std::string &outName,
                                      double &thickness);

  /** @brief periodicMeshMapper method generates map of nodes on periodic mesh
             interfaces at boundaries for FEM analysis.
      @param patch1 First boundary for mapping.
      @param patch2 Second boundary for mapping.
  **/
  void periodicMeshMapper(std::string &patch1, std::string &patch2);

  // --- private methods
 private:
  /** @brief createPatchDict utility creates a dictionary needed for createPatch
            method.
      @param dirStat domain number to skip (output from splitMshRegions)
      @param write Write boolean
  **/
  void createPatchDict(const int &dirStat, const bool &write);

  // Pointer for parameters
 private:
  /** @brief MeshManipulationFoam object
   **/
  MeshManipulationFoamParams *_mshMnipPrms;
  std::unique_ptr<Foam::dictionary> cpDict_;
};

#endif  // NEMOSYS_MESHMANIPULATIONFOAM_H_
