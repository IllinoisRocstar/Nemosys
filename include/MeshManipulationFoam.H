#ifndef NEMOSYS_MESHMANIPULATIONFOAM_H_
#define NEMOSYS_MESHMANIPULATIONFOAM_H_

#include <vtkDoubleArray.h>
#include <vtkMeshQuality.h>

#include "nemosys_export.h"

// OpenFOAM
#include "argList.H"
#include "fileName.H"
#include "fvCFD.H"
#include "fvMesh.H"

// SurfLambdaMuSmooth
#include "MeshedSurfaces.H"
#include "argList.H"
#include "boundBox.H"
#include "edgeMesh.H"
#include "matchPoints.H"

// splitMeshByRegions
#include "EdgeMap.H"
#include "IOobjectList.H"
#include "ReadFields.H"
#include "SortableList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "fvMeshSubset.H"
#include "fvMeshTools.H"
#include "mappedWallPolyPatch.H"
#include "polyTopoChange.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "syncTools.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// mergeMeshes
#include "Time.H"
#include "mergePolyMesh.H"

// createPatch
#include "IOPtrList.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "SortableList.H"
#include "cyclicPolyPatch.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "polyModifyFace.H"
#include "syncTools.H"
#include "wordReList.H"

// Nemosys
#include "MeshManipulationFoamParams.H"
#include "meshBase.H"

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
  MeshManipulationFoam(MeshManipulationFoamParams* params) {
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
  std::pair<int, int> splitMshRegions();  // SplitMeshRegions utility

  /** @brief .
      @param dirStat domain number to skip (output from splitMshRegions)
      @param nDomains Number of meshes to merge
  **/
  void mergeMeshes(int dirStat, int nDomains);  // MergeMeshes utility

  /** @brief createPatch utility creates user-defined patches in Foam mesh. User
            input is usually createPatchDict.
      @param dirStat domain number to skip (output from splitMshRegions)
  **/
  void createPatch(int dirStat);  // CreatePatch utility

  /** @brief foamToSurface utility reads Foam mesh, extracts its surface and
            writes into an STL file.
  **/
  void foamToSurface();  // FoamToSurface utility

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
  void addCohesiveElements(double tol, const std::string outName);

  /** @brief createControlDict method creates a general purpose controlDict for
            Foam environment.
  **/
  void createControlDict();

  /** @brief This method adds connectivity elements with finite thickness
        at shared interface between two conformal meshes
      @param tol Tolerance for conformal duplicate nodes between two meshes.
      @param outName Output mesh name with extension (i.e "output.vtu")
      @param thickness Desired thickness for new cohesive elements.
  **/
  void addArtificialThicknessElements(double tol, const std::string outName,
                                      double thickness);

  /** @brief periodicMeshMapper method generates map of nodes on periodic mesh
             interfaces at boundaries for FEM analysis.
      @param patch1 First boundary for mapping.
      @param patch2 Second boundary for mapping.
  **/
  void periodicMeshMapper(std::string patch1, std::string patch2);

  /** @brief Initializes AMR workflow if called. (Used in RefineDriver)
      @param startT Simulation start time.
      @param dt Time step for simulation
      @param endT End time for simulation
  **/
  void initAMRWorkFlow(const double startT, const double& dt,
                       const double& endT);

  // --- private methods
 private:
  /** @brief createPatchDict utility creates a dictionary needed for createPatch
            method.
      @param dirStat domain number to skip (output from splitMshRegions)
  **/
  void createPatchDict(int dirStat);  // Creates createPatchDict

  /** @brief createVtkCell method created a cell in VTK database.
      @param dataSet VTK mesh database
      @param cellType VTK cell type number
      @param vrtIds Cell bounding vertices
  **/
  void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
                     const int cellType, std::vector<int>& vrtIds);

  // surfaceLambdaMuSmooth functions
  /** @brief surfLambdaMuSmooth Internal Function. Calculates averages
            between fixed points
  **/
  Foam::tmp<Foam::pointField> avg(const Foam::meshedSurface& s,
                                  const Foam::PackedBoolList& fixedPoints);

  /** @brief surfLambdaMuSmooth Internal Function. Obtains fixed mesh points
            from edges
  **/
  void getFixedPoints(const Foam::edgeMesh& feMesh,
                      const Foam::pointField& points,
                      Foam::PackedBoolList& fixedPoints);

  // splitMshByRegion functions
  /** @brief splitMshByRegion Internal Function. Renames patches discovered
            by walking in mesh
  **/
  void renamePatches(Foam::fvMesh& mesh, const Foam::word& prefix,
                     const Foam::labelList& patchesToRename);

  /** @brief splitMshByRegion Internal Function. Makes volVectorFields of
            subset volumes
  **/
  template <class GeoField>
  void subsetVolFields(const Foam::fvMesh& mesh, const Foam::fvMesh& subMesh,
                       const Foam::labelList& cellMap,
                       const Foam::labelList& faceMap,
                       const Foam::labelHashSet& addedPatches);

  /** @brief splitMshByRegion Internal Function. Makes volSurfaceFields of
            subset surfaces
  **/
  template <class GeoField>
  void subsetSurfaceFields(const Foam::fvMesh& mesh,
                           const Foam::fvMesh& subMesh,
                           const Foam::labelList& cellMap,
                           const Foam::labelList& faceMap,
                           const Foam::labelHashSet& addedPatches);

  /** @brief splitMshByRegion Internal Function. Outputs non-region
            disconnected cells
  **/
  Foam::labelList getNonRegionCells(const Foam::labelList& cellRegion,
                                    const Foam::label regionI);

  /** @brief splitMshByRegion Internal Function. Adds new region to interface.
   **/
  void addToInterface(const Foam::polyMesh& mesh, const Foam::label zoneID,
                      const Foam::label ownRegion, const Foam::label neiRegion,
                      Foam::EdgeMap<Foam::Map<Foam::label>>& regionsToSize);

  /** @brief splitMshByRegion Internal Function. Interface size in terms of
            cells
  **/
  void getInterfaceSizes(const Foam::polyMesh& mesh, const bool useFaceZones,
                         const Foam::labelList& cellRegion,
                         const Foam::wordList& regionNames,
                         Foam::edgeList& interfaces,
                         Foam::List<Foam::Pair<Foam::word>>& interfaceNames,
                         Foam::labelList& interfaceSizes,
                         Foam::labelList& faceToInterface);

  /** @brief splitMshByRegion Internal Function. Defines mesh of disconnected
            regions after separating
  **/
  Foam::autoPtr<Foam::mapPolyMesh> createRegionMesh(
      const Foam::fvMesh& mesh, const Foam::labelList& cellRegion,
      const Foam::label regionI, const Foam::word& regionName,
      const Foam::labelList& interfacePatches,
      const Foam::labelList& faceToInterface,
      Foam::autoPtr<Foam::fvMesh>& newMesh);

  /** @brief splitMshByRegion Internal Function. Writes region mesh into
            separate folders
  **/
  void createAndWriteRegion(const Foam::fvMesh& mesh,
                            const Foam::labelList& cellRegion,
                            const Foam::wordList& regionNames,
                            const bool prefixRegion,
                            const Foam::labelList& faceToInterface,
                            const Foam::labelList& interfacePatches,
                            const Foam::label regionI,
                            const Foam::word& newMeshInstance);

  /** @brief splitMshByRegion Internal Function. Adds patches for new
            regions to mesh
  **/
  Foam::labelList addRegionPatches(
      Foam::fvMesh& mesh, const Foam::wordList& regionNames,
      const Foam::edgeList& interfaces,
      const Foam::List<Foam::Pair<Foam::word>>& interfaceNames);

  /** @brief splitMshByRegion Internal Function. Finds certain region for
            writing
  **/
  Foam::label findCorrespondingRegion(const Foam::labelList& existingZoneID,
                                      const Foam::labelList& cellRegion,
                                      const Foam::label nCellRegions,
                                      const Foam::label zoneI,
                                      const Foam::label minOverlapSize);

  /** @brief splitMshByRegion Internal Function. Outputs zoneID
   **/
  void getZoneID(const Foam::polyMesh& mesh,
                 const Foam::cellZoneMesh& cellZones, Foam::labelList& zoneID,
                 Foam::labelList& neiZoneID);

  /** @brief splitMshByRegion Internal Function. Method for sloppyCellZone
            features
  **/
  int matchRegions(const bool sloppyCellZones, const Foam::polyMesh& mesh,
                   const Foam::label nCellRegions,
                   const Foam::labelList& cellRegion,
                   Foam::labelList& regionToZone, Foam::wordList& regionNames,
                   Foam::labelList& zoneToRegion);

  /** @brief splitMshByRegion Internal Function. PolyMesh writing handler
   **/
  void writeCellToRegion(const Foam::fvMesh& mesh,
                         const Foam::labelList& cellRegion);

  // mergeMeshes functions
  /** @brief mergeMeshes Internal Function
   **/
  void getRootCase(Foam::fileName& casePath);

  // createPatch functions
  /** @brief createPatch Internal Function. Changes patch id's in mesh.
   **/
  void changePatchID(const Foam::polyMesh& mesh, const Foam::label faceID,
                     const Foam::label patchID, Foam::polyTopoChange& meshMod);

  /** @brief createPatch Internal Function. Looks for duplicate and
            non-zero patches
  **/
  void filterPatches(Foam::polyMesh& mesh,
                     const Foam::HashSet<Foam::word>& addedPatchNames);

  /** @brief createPatch Internal Function. Moving mesh function
   **/
  void dumpCyclicMatch(const Foam::fileName& prefix,
                       const Foam::polyMesh& mesh);

  /** @brief createPatch Internal Function. Lists all patches to be merged
   **/
  void separateList(const Foam::vectorField& separation,
                    Foam::UList<Foam::vector>& field);

  /** @brief createPatch Internal Function. Ensures conformality
            if pointSync option is true
  **/
  template <class CombineOp>
  void syncPoints(const Foam::polyMesh& mesh, Foam::pointField& points,
                  const CombineOp& cop, const Foam::point& nullValue);

  // Pointer for parameters
 private:
  /** @brief MeshManipulationFoam object
   **/
  MeshManipulationFoamParams* _mshMnipPrms;
};

#endif  // NEMOSYS_MESHMANIPULATIONFOAM_H_
