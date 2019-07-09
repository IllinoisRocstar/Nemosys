#ifndef MESHMANIPULATIONFOAM_H
#define MESHMANIPULATIONFOAM_H

#include <vtkMeshQuality.h>
#include <vtkDoubleArray.h>

#include "MeshManipulationFoamParams.H"
#include "meshBase.H"
#include "argList.H"
#include "polyMeshGenModifier.H"
#include "meshOptimizer.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "fileName.H"

// SurfLambdaMuSmooth
#include "argList.H"
#include "boundBox.H"
#include "edgeMesh.H"
#include "matchPoints.H"
#include "MeshedSurfaces.H"

// splitMeshByRegions
#include "SortableList.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "EdgeMap.H"
#include "syncTools.H"
#include "ReadFields.H"
#include "mappedWallPolyPatch.H"
#include "fvMeshTools.H"
#include "zeroGradientFvPatchFields.H"

// mergeMeshes
#include "Time.H"
#include "mergePolyMesh.H"

// createPatch
#include "cyclicPolyPatch.H"
#include "syncTools.H"
#include "polyMesh.H"
#include "SortableList.H"
#include "OFstream.H"
#include "meshTools.H"
#include "IOPtrList.H"
#include "polyModifyFace.H"
#include "wordReList.H"
#include "IOdictionary.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

/*
  MeshManipulation class contains several OpenFOAM utilities to perform
  various mesh manipulation operations. These utilities reads Foam mesh
  and outputs Foam mesh. Some of them are surface mesh manipulation
  utilities.
*/
class MeshManipulationFoamParams; // Friend Class

class MeshManipulationFoam
{
  public:
  MeshManipulationFoam(){}; // Standard Constructor
  MeshManipulationFoam(MeshManipulationFoamParams* params)
  {_mshMnipPrms = params;}  // Alternate Constructor
  ~MeshManipulationFoam();  // Class Destructor


  public:
  void createControlDict();
  void surfLambdaMuSmooth();  // SurfaceLambdaMuSmooth utility
  int splitMshRegions();    // SplitMeshRegions utility
  void mergeMeshes(int dirStat, int nDomains);  // MergeMeshes utility
  void createPatch();   // CreatePatch utility
  void foamToSurface();   // FoamToSurface utility
  int surfSpltByTopology();  // surfaceSplitByTopology utility
  //void addConnectivityData(); // Adds connectivity info at duplicate nodes.


  private:
  void createPatchDict();   // Creates createPatchDict
  //void createVtkCell(vtkSmartPointer<vtkUnstructuredGrid> dataSet,
  //      const int cellType, std::vector<int>& vrtIds);

  // surfaceLambdaMuSmooth functions
  // Defines average between points
  Foam::tmp<Foam::pointField> avg(const Foam::meshedSurface& s,
    const Foam::PackedBoolList& fixedPoints);

  // Obtaines fixed mesh points from edges  
  void getFixedPoints(const Foam::edgeMesh& feMesh,
    const Foam::pointField& points,
    Foam::PackedBoolList& fixedPoints);

  // splitMeshByRegion functions
  // Renames patches discovered by walking in mesh
  void renamePatches(Foam::fvMesh& mesh,
    const Foam::word& prefix,
    const Foam::labelList& patchesToRename);

  // Makes volVectorFields of subset volumes
  template<class GeoField>
  void subsetVolFields(const Foam::fvMesh& mesh,
    const Foam::fvMesh& subMesh,
    const Foam::labelList& cellMap,
    const Foam::labelList& faceMap,
    const Foam::labelHashSet& addedPatches);

  // Makes volSurfaceFields of subset surfaces
  template<class GeoField>
  void subsetSurfaceFields(const Foam::fvMesh& mesh,
    const Foam::fvMesh& subMesh,
    const Foam::labelList& cellMap,
    const Foam::labelList& faceMap,
    const Foam::labelHashSet& addedPatches);

  // Outputs non-region disconnected cells
  Foam::labelList getNonRegionCells(const Foam::labelList& cellRegion,
    const Foam::label regionI);

  // Adds new region to interface
  void addToInterface(const Foam::polyMesh& mesh,
    const Foam::label zoneID, const Foam::label ownRegion,
    const Foam::label neiRegion,
    Foam::EdgeMap<Foam::Map<Foam::label>>& regionsToSize);

  // Interface size in terms of cells
  void getInterfaceSizes(const Foam::polyMesh& mesh,
    const bool useFaceZones, const Foam::labelList& cellRegion,
    const Foam::wordList& regionNames,
    Foam::edgeList& interfaces,
    Foam::List<Foam::Pair<Foam::word>>& interfaceNames,
    Foam::labelList& interfaceSizes,
    Foam::labelList& faceToInterface);

  // Defines mesh of disconnected regions after separating
  Foam::autoPtr<Foam::mapPolyMesh> createRegionMesh(
    const Foam::fvMesh& mesh, const Foam::labelList& cellRegion,
    const Foam::label regionI, const Foam::word& regionName,
    const Foam::labelList& interfacePatches,
    const Foam::labelList& faceToInterface,
    Foam::autoPtr<Foam::fvMesh>& newMesh);
  
  // Writes region mesh into separate folders
  void createAndWriteRegion(const Foam::fvMesh& mesh,
    const Foam::labelList& cellRegion,
    const Foam::wordList& regionNames,
    const bool prefixRegion,
    const Foam::labelList& faceToInterface,
    const Foam::labelList& interfacePatches,
    const Foam::label regionI,
    const Foam::word& newMeshInstance);

  // Adds patches for new regions to mesh
  Foam::labelList addRegionPatches(Foam::fvMesh& mesh,
    const Foam::wordList& regionNames,
    const Foam::edgeList& interfaces,
    const Foam::List<Foam::Pair<Foam::word>>& interfaceNames);

  // Finds certain region for writing
  Foam::label findCorrespondingRegion(
    const Foam::labelList& existingZoneID,
    const Foam::labelList& cellRegion,
    const Foam::label nCellRegions,
    const Foam::label zoneI,
    const Foam::label minOverlapSize);
  
  // Outputs zoneID
  void getZoneID(
    const Foam::polyMesh& mesh,
    const Foam::cellZoneMesh& cellZones,
    Foam::labelList& zoneID,
    Foam::labelList& neiZoneID);

  // Method for sloppyCellZone features
  int matchRegions(
    const bool sloppyCellZones,
    const Foam::polyMesh& mesh,
    const Foam::label nCellRegions,
    const Foam::labelList& cellRegion,
    Foam::labelList& regionToZone,
    Foam::wordList& regionNames,
    Foam::labelList& zoneToRegion);

  // PolyMesh writing handler
  void writeCellToRegion(
    const Foam::fvMesh& mesh, const Foam::labelList& cellRegion);

  //mergeMeshes functions
  void getRootCase(Foam::fileName& casePath);

  //createPatch functions
  void changePatchID(const Foam::polyMesh& mesh,
    const Foam::label faceID, const Foam::label patchID,
    Foam::polyTopoChange& meshMod);

  // Looks for duplicate and non-zero patches
  void filterPatches(Foam::polyMesh& mesh,
    const Foam::HashSet<Foam::word>& addedPatchNames);

  // Moving mesh function
  void dumpCyclicMatch(
    const Foam::fileName& prefix, 
    const Foam::polyMesh& mesh);

  // Lists all patches to be merged
  void separateList(
    const Foam::vectorField& separation,
    Foam::UList<Foam::vector>& field);

  // Ensures conformality if pointSync option is true
  template<class CombineOp>
  void syncPoints(const Foam::polyMesh& mesh,
    Foam::pointField& points,
    const CombineOp& cop,
    const Foam::point& nullValue);

  // Pointer for parameters
  private:
  MeshManipulationFoamParams* _mshMnipPrms;
  //vtkSmartPointer<vtkUnstructuredGrid> packDS;
  //vtkSmartPointer<vtkUnstructuredGrid> surroundingDS;
};

#endif