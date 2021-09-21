#ifndef SNAPPYMESH_H
#define SNAPPYMESH_H

#include "IOmanip.H"
#include "MeshedSurface.H"
#include "Time.H"
#include "UnsortedMeshedSurface.H"
#include "argList.H"
#include "cellModeller.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "faceSet.H"
#include "fileName.H"
#include "fvMesh.H"
#include "fvMeshDistribute.H"
#include "fvMeshTools.H"
#include "getDicts.H"
#include "globalIndex.H"
#include "layerParameters.H"
#include "meshRefinement.H"
#include "motionSmoother.H"
#include "noDecomp.H"
#include "polyTopoChange.H"
#include "processorMeshes.H"
#include "profiling.H"
#include "refinementFeatures.H"
#include "refinementParameters.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "shellSurfaces.H"
#include "snapParameters.H"
#include "snappyLayerDriver.H"
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyVoxelMeshDriver.H"
#include "surfZoneIdentifierList.H"
#include "uindirectPrimitivePatch.H"
#include "vtkSetWriter.H"
#include "wallPolyPatch.H"

class snappyMesh {
 public:
  /** @brief Constructor */
  snappyMesh();

  /** @brief Destructor */
  ~snappyMesh();

  /**
   * @brief Generate new mesh using snappyHexMeshDict (in file directory or in
   * memory)
   */
  Foam::fvMesh *generate(
      const std::unique_ptr<Foam::Time> &runTime = NULL,
      const std::unique_ptr<Foam::dictionary> &snappyMshDict_ = NULL,
      const std::unique_ptr<Foam::fvMesh> &fmesh_ = NULL,
      const bool &write = true);

 private:
  /**
   * @brief Convert size (as fraction of defaultCellSize) to refinement level
   * @param level0Coeff ratio of hex cell size v.s. defaultCellSize
   * @param sizeCoeff
   */
  Foam::label sizeCoeffToRefinement(const Foam::scalar level0Coeff,
                                    const Foam::scalar sizeCoeff);

  Foam::autoPtr<Foam::refinementSurfaces> createRefinementSurfaces(
      const Foam::searchableSurfaces &allGeometry,
      const Foam::dictionary &surfacesDict,
      const Foam::dictionary &shapeControlDict,
      const Foam::label gapLevelIncrement, const Foam::scalar level0Coeff);

  void extractSurface(const Foam::polyMesh &mesh, const Foam::Time &runTime,
                      const Foam::labelHashSet &includePatches,
                      const Foam::fileName &outFileName);

  Foam::label checkAlignment(const Foam::polyMesh &mesh, const Foam::scalar tol,
                             Foam::Ostream &os);

  /** @brief Check writing tolerance before doing any serious work */
  Foam::scalar getMergeDistance(const Foam::polyMesh &mesh,
                                const Foam::scalar mergeTol, const bool dryRun);

  void removeZeroSizedPatches(Foam::fvMesh &mesh);

  /** Write mesh and additional information */
  void writeMesh(const std::string &msg,
                 const Foam::meshRefinement &meshRefiner,
                 const Foam::meshRefinement::debugType debugLevel,
                 const Foam::meshRefinement::writeType writeLevel);
};

#endif  // SNAPPYMESH_H
