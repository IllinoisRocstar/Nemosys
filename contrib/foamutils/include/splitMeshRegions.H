#ifndef SPLITMESHREGIONS_H
#define SPLITMESHREGIONS_H

#include "EdgeMap.H"
#include "IOobjectList.H"
#include "ReadFields.H"
#include "SortableList.H"
#include "Time.H"
#include "argList.H"
#include "cellSet.H"
#include "faceSet.H"
#include "fileName.H"
#include "fvMesh.H"
#include "fvMeshSubset.H"
#include "fvMeshTools.H"
#include "mappedWallPolyPatch.H"
#include "polyTopoChange.H"
#include "regionSplit.H"
#include "removeCells.H"
#include "syncTools.H"
#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

class splitMeshRegions {
 public:
  /** @brief Constructor */
  splitMeshRegions();

  /** @brief Destructor */
  ~splitMeshRegions();

  /**
   * @brief Splits the mesh based on cellZones or regions
   * @param useCellZones Use cell zones boolean
   * @return pair of vectors with some information regarding packmesh
   */
  std::pair<std::vector<int>, std::vector<std::string>> execute(
      const bool useCellZones);

 private:
  void renamePatches(Foam::fvMesh &mesh, const Foam::word &prefix,
                     const Foam::labelList &patchesToRename);

  template <class GeoField>
  void subsetVolFields(const Foam::fvMesh &mesh, const Foam::fvMesh &subMesh,
                       const Foam::labelList &cellMap,
                       const Foam::labelList &faceMap,
                       const Foam::labelHashSet &addedPatches);

  template <class GeoField>
  void subsetSurfaceFields(const Foam::fvMesh &mesh,
                           const Foam::fvMesh &subMesh,
                           const Foam::labelList &cellMap,
                           const Foam::labelList &faceMap,
                           const Foam::labelHashSet &addedPatches);

  Foam::labelList getNonRegionCells(const Foam::labelList &cellRegion,
                                    const Foam::label regionI);

  void addToInterface(const Foam::polyMesh &mesh, const Foam::label zoneID,
                      const Foam::label ownRegion, const Foam::label neiRegion,
                      Foam::EdgeMap<Foam::Map<Foam::label>> &regionsToSize);

  void getInterfaceSizes(const Foam::polyMesh &mesh, const bool useFaceZones,
                         const Foam::labelList &cellRegion,
                         const Foam::wordList &regionNames,
                         Foam::edgeList &interfaces,
                         Foam::List<Foam::Pair<Foam::word>> &interfaceNames,
                         Foam::labelList &interfaceSizes,
                         Foam::labelList &faceToInterface);

  Foam::autoPtr<Foam::mapPolyMesh> createRegionMesh(
      const Foam::fvMesh &mesh, const Foam::labelList &cellRegion,
      const Foam::label regionI, const Foam::word &regionName,
      const Foam::labelList &interfacePatches,
      const Foam::labelList &faceToInterface,
      Foam::autoPtr<Foam::fvMesh> &newMesh);

  void createAndWriteRegion(const Foam::fvMesh &mesh,
                            const Foam::labelList &cellRegion,
                            const Foam::wordList &regionNames,
                            const bool prefixRegion,
                            const Foam::labelList &faceToInterface,
                            const Foam::labelList &interfacePatches,
                            const Foam::label regionI,
                            const Foam::word &newMeshInstance);

  Foam::labelList addRegionPatches(
      Foam::fvMesh &mesh, const Foam::wordList &regionNames,
      const Foam::edgeList &interfaces,
      const Foam::List<Foam::Pair<Foam::word>> &interfaceNames);

  Foam::label findCorrespondingRegion(
      const Foam::labelList &existingZoneID,  // per cell the (unique) zoneID
      const Foam::labelList &cellRegion, const Foam::label nCellRegions,
      const Foam::label zoneI, const Foam::label minOverlapSize);

  void getZoneID(const Foam::polyMesh &mesh,
                 const Foam::cellZoneMesh &cellZones, Foam::labelList &zoneID,
                 Foam::labelList &neiZoneID);

  int matchRegions(const bool sloppyCellZones, const Foam::polyMesh &mesh,

                   const Foam::label nCellRegions,
                   const Foam::labelList &cellRegion,

                   Foam::labelList &regionToZone, Foam::wordList &regionNames,
                   Foam::labelList &zoneToRegion);

  void writeCellToRegion(const Foam::fvMesh &mesh,
                         const Foam::labelList &cellRegion);
};

#endif  // SPLITMESHREGIONS_H
