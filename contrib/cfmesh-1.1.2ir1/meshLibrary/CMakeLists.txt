# Finding libraries ############################################################
find_package(OpenMP REQUIRED)

# Setting library files ########################################################

# Setting include files
file(GLOB_RECURSE CFMESH_INC_FILES *.H)

set(CFMESH_INC_DIRS "")
foreach(file_path ${CFMESH_INC_FILES})
  get_filename_component(dir_path ${file_path} PATH)
  list(APPEND CFMESH_INC_DIRS ${dir_path})
endforeach()
list(REMOVE_DUPLICATES CFMESH_INC_DIRS)

# Setting cfMesh source
set(CFMESH_SRCS
utilities/checkMeshDict/checkMeshDict.C

utilities/containers/Lists/pointFieldPMG.C
utilities/containers/Lists/faceListPMG.C

utilities/containers/VRWGraph/VRWGraph.C
utilities/containers/VRWGraph/VRWGraphSMPModifier.C

utilities/containers/VRWGraphList/VRWGraphList.C

utilities/containers/Graphs/cellIOGraph.C
utilities/containers/Graphs/faceIOGraph.C

utilities/meshes/polyMeshGen/polyMeshGen.C
utilities/meshes/polyMeshGen/polyMeshGenPoints.C
utilities/meshes/polyMeshGen/polyMeshGenFaces.C
utilities/meshes/polyMeshGen/polyMeshGenCells.C

utilities/meshes/polyMeshGen2DEngine/polyMeshGen2DEngine.C

utilities/anisotropicMeshing/coordinateModification/coordinateModification.C
utilities/anisotropicMeshing/coordinateModification/newCoordinateModification.C
utilities/anisotropicMeshing/coordinateModification/planeScaling.C
utilities/anisotropicMeshing/coordinateModification/boxScaling.C

utilities/anisotropicMeshing/edgeMeshGeometryModification/edgeMeshGeometryModification.C
utilities/anisotropicMeshing/surfaceMeshGeometryModification/surfaceMeshGeometryModification.C
utilities/anisotropicMeshing/polyMeshGenGeometryModification/polyMeshGenGeometryModification.C

utilities/anisotropicMeshing/coordinateModification/coordinateModifier.C

utilities/meshes/polyMeshGen/boundaryPatch/boundaryPatchBase.C
utilities/meshes/polyMeshGen/boundaryPatch/boundaryPatch.C
utilities/meshes/polyMeshGen/boundaryPatch/processorBoundaryPatch.C

utilities/meshes/polyMeshGenModifier/polyMeshGenModifierRemoveUnusedVertices.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierRemoveFaces.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierRemoveCells.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierReorderBoundaryFaces.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierAddCells.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierAddProcessorFaces.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierAddBufferCells.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierReplaceBoundary.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierZipUpCells.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierRenumberMesh.C
utilities/meshes/polyMeshGenModifier/polyMeshGenModifierAddCellByCell.C

utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingCellCells.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingCellEdges.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingCellPoints.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingCentresAndAreas.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingCentresAndVols.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingClear.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingUpdateGeometry.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingEdgeCells.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingEdgeFaces.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingFaceEdges.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingEdges.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingPointCells.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingPointEdges.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingPointFaces.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingPointPoints.C
utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressingParallelAddressing.C

utilities/meshes/polyMeshGenChecks/polyMeshGenChecks.C
utilities/meshes/polyMeshGenChecks/polyMeshGenChecksGeometry.C
utilities/meshes/polyMeshGenChecks/polyMeshGenChecksTopology.C

utilities/meshes/partTetMesh/partTetMesh.C
utilities/meshes/partTetMesh/partTetMeshAddressing.C
utilities/meshes/partTetMesh/partTetMeshParallelAddressing.C
utilities/meshes/partTetMesh/partTetMeshSimplex.C

utilities/meshes/partTriMesh/partTriMesh.C
utilities/meshes/partTriMesh/partTriMeshAddressing.C
utilities/meshes/partTriMesh/partTriMeshParallelAddressing.C
utilities/meshes/partTriMesh/partTriMeshSimplex.C

utilities/meshes/triSurf/triSurf.C
utilities/meshes/triSurf/triSurfPoints.C
utilities/meshes/triSurf/triSurfFacets.C
utilities/meshes/triSurf/triSurfFeatureEdges.C
utilities/meshes/triSurf/triSurfAddressing.C
utilities/meshes/triSurf/triSurfModifier.C

utilities/intersectionTools/findCellsIntersectingSurface/findCellsIntersectingSurface.C

utilities/tetrahedra/tetCreatorOctree/tetCreatorOctree.C
utilities/tetrahedra/tetCreatorOctree/tetCreatorOctreePointsAndAddressing.C
utilities/tetrahedra/tetCreatorOctree/tetCreatorOctreeFromFacesWithCentreNode.C
utilities/tetrahedra/tetCreatorOctree/tetCreatorOctreeTetsAroundEdges.C
utilities/tetrahedra/tetCreatorOctree/tetCreatorOctreeTetsAroundSplitEdges.C
utilities/tetrahedra/tetCreatorOctree/tetCreatorOctreeTetsFromSplitFaces.C

utilities/faceDecomposition/faceDecomposition.C
utilities/faceDecomposition/decomposeFaces.C

utilities/helperFunctions/helperFunctionsStringConversion.C

utilities/helperClasses/sortEdgesIntoChains/sortEdgesIntoChains.C

utilities/surfaceTools/surfaceMorpherCells/surfaceMorpherCells.C
utilities/surfaceTools/surfaceMorpherCells/surfaceMorpherCellsMorphInternalFaces.C
utilities/surfaceTools/surfaceMorpherCells/surfaceMorpherCellsCreateBoundaryFaces.C

utilities/decomposeCells/decomposeCells.C
utilities/decomposeCells/decomposeCellsPyramids.C
utilities/decomposeCells/decomposeCellsDecomposition.C

utilities/smoothers/topology/topologicalCleaner/topologicalCleaner.C
utilities/smoothers/topology/topologicalCleaner/topologyCleanerNonConsecutiveBoundaryVertices.C
utilities/smoothers/topology/topologicalCleaner/topologicalCleanerInvalidVertices.C
utilities/smoothers/topology/topologicalCleaner/topologyCleanerNonMappableCells.C
utilities/smoothers/topology/checkCellConnectionsOverFaces/checkCellConnectionsOverFaces.C
utilities/smoothers/topology/checkIrregularSurfaceConnections/checkIrregularSurfaceConnections.C
utilities/smoothers/topology/checkIrregularSurfaceConnections/checkIrregularSurfaceConnectionsFunctions.C
utilities/smoothers/topology/checkNonMappableCellConnections/checkNonMappableCellConnections.C
utilities/smoothers/topology/checkBoundaryFacesSharingTwoEdges/checkBoundaryFacesSharingTwoEdges.C

utilities/boundaryLayers/boundaryLayers/boundaryLayers.C
utilities/boundaryLayers/boundaryLayers/boundaryLayersCreateVertices.C
utilities/boundaryLayers/boundaryLayers/boundaryLayersFacesAndCells.C
utilities/boundaryLayers/boundaryLayers/boundaryLayerCells.C
utilities/boundaryLayers/boundaryLayers/boundaryLayersCheckTopologyOfBndFaces.C
utilities/boundaryLayers/boundaryLayers/boundaryLayersWrapperLayer.C
utilities/boundaryLayers/extrudeLayer/extrudeLayer.C
utilities/boundaryLayers/detectBoundaryLayers/detectBoundaryLayers.C
utilities/boundaryLayers/detectBoundaryLayers/detectBoundaryLayersFunctions.C
utilities/boundaryLayers/refineBoundaryLayers/refineBoundaryLayers.C
utilities/boundaryLayers/refineBoundaryLayers/refineBoundaryLayersFunctions.C
utilities/boundaryLayers/refineBoundaryLayers/refineBoundaryLayersFaces.C
utilities/boundaryLayers/refineBoundaryLayers/refineBoundaryLayersCells.C
utilities/boundaryLayers/triangulateNonPlanarBaseFaces/triangulateNonPlanarBaseFaces.C
utilities/boundaryLayers/triangulateNonPlanarBaseFaces/triangulateNonPlanarBaseFacesFunctions.C

utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.C
utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngineCalculateBoundaryNodesAndFaces.C
utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngineParallelAddressing.C
utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngineModifier.C

utilities/surfaceTools/meshSurfacePartitioner/meshSurfacePartitioner.C
utilities/surfaceTools/meshSurfacePartitioner/meshSurfacePartitionerFunctions.C

utilities/surfaceTools/meshSurfaceCheckInvertedVertices/meshSurfaceCheckInvertedVertices.C

utilities/surfaceTools/meshSurfaceCheckEdgeTypes/meshSurfaceCheckEdgeTypes.C

utilities/surfaceTools/meshSurfaceMapper/meshSurfaceMapper.C
utilities/surfaceTools/meshSurfaceMapper/meshSurfaceMapperMapVertices.C
utilities/surfaceTools/meshSurfaceMapper/meshSurfaceMapperCornersAndEdges.C
utilities/surfaceTools/meshSurfaceMapper/meshSurfaceMapperPremapVertices.C

utilities/surfaceTools/meshSurfaceMapper2D/meshSurfaceMapper2D.C
utilities/surfaceTools/meshSurfaceMapper2D/meshSurfaceMapper2DMapVertices.C
utilities/surfaceTools/meshSurfaceMapper2D/meshSurfaceMapper2DPremapVertices.C

utilities/surfaceTools/edgeExtraction/edgeExtractor/edgeExtractor.C
utilities/surfaceTools/edgeExtraction/edgeExtractor/edgeExtractorCorners.C

utilities/surfaceTools/meshSurfaceEdgeExtractorNonTopo/meshSurfaceEdgeExtractorNonTopo.C
utilities/surfaceTools/meshSurfaceEdgeExtractorNonTopo/meshSurfaceEdgeExtractorNonTopoDistributeFaces.C

utilities/surfaceTools/meshSurfaceEdgeExtractor2D/meshSurfaceEdgeExtractor2D.C
utilities/surfaceTools/meshSurfaceEdgeExtractor2D/meshSurfaceEdgeExtractor2DDistributeFaces.C

utilities/surfaceTools/meshSurfaceEdgeExtractorFUN/meshSurfaceEdgeExtractorFUN.C
utilities/surfaceTools/meshSurfaceEdgeExtractorFUN/meshSurfaceEdgeExtractorFUNDistributeFaces.C

utilities/surfaceTools/createFundamentalSheets/createFundamentalSheets/createFundamentalSheets.C
utilities/surfaceTools/createFundamentalSheets/createFundamentalSheetsFJ/createFundamentalSheetsFJ.C
utilities/surfaceTools/createFundamentalSheets/createFundamentalSheetsJFS/createFundamentalSheetsJFS.C

utilities/surfaceTools/correctEdgesBetweenPatches/correctEdgesBetweenPatches.C
utilities/surfaceTools/correctEdgesBetweenPatches/correctEdgesBetweenPatchesDistributeFaces.C

utilities/surfaceTools/renameBoundaryPatches/renameBoundaryPatches.C

utilities/smoothers/geometry/meshOptimizer/meshOptimizer.C
utilities/smoothers/geometry/meshOptimizer/meshOptimizerOptimizePoint.C
utilities/smoothers/geometry/meshOptimizer/meshOptimizerOptimizePointParallel.C
utilities/smoothers/geometry/meshOptimizer/meshOptimizerOptimizeSurface.C
utilities/smoothers/geometry/meshOptimizer/optimizeMeshFV.C

utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/tetMeshOptimisation.C
utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/tetMeshOptimisationParallel.C

utilities/smoothers/geometry/meshOptimizer/symmetryPlaneOptimisation/symmetryPlaneOptimisation.C

utilities/smoothers/geometry/meshOptimizer/boundaryLayerOptimisation/boundaryLayerOptimisation.C
utilities/smoothers/geometry/meshOptimizer/boundaryLayerOptimisation/boundaryLayerOptimisationFunctions.C
utilities/smoothers/geometry/meshOptimizer/boundaryLayerOptimisation/boundaryLayerOptimisationNormals.C
utilities/smoothers/geometry/meshOptimizer/boundaryLayerOptimisation/boundaryLayerOptimisationThickness.C

utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/simplexSmoother/simplexSmoother.C

utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/knuppMetric/knuppMetric.C

utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/meshUntangler/meshUntangler.C
utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/meshUntangler/meshUntanglerCutRegion.C
utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/meshUntangler/meshUntanglerCutRegionPoints.C
utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/meshUntangler/meshUntanglerCutRegionEdges.C
utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/meshUntangler/meshUntanglerCutRegionFaces.C
utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/meshUntangler/meshUntanglerCutRegionTieBreak.C

utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/quadricMetric/quadricMetric.C

utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/volumeOptimizer/volumeOptimizer.C
utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/volumeOptimizer/volumeOptimizerEvaluateGradients.C

utilities/smoothers/geometry/meshSurfaceOptimizer/meshSurfaceOptimizer.C
utilities/smoothers/geometry/meshSurfaceOptimizer/meshSurfaceOptimizerCalculateTrianglesAndAddressing.C
utilities/smoothers/geometry/meshSurfaceOptimizer/meshSurfaceOptimizerOptimizePoint.C
utilities/smoothers/geometry/meshSurfaceOptimizer/meshSurfaceOptimizerOptimizeSurface.C
utilities/smoothers/geometry/meshSurfaceOptimizer/meshSurfaceOptimizerOptimizePointParallel.C

utilities/smoothers/geometry/meshSurfaceOptimizer/advancedSurfaceSmoothers/surfaceOptimizer/surfaceOptimizer.C

utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCube.C
utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeIntersections.C
utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeRecursiveFunctions.C
utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeRefine.C
utilities/octrees/meshOctree/meshOctreeCube/meshOctreeCubeCoordinatesIntersections.C

utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifier.C
utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierRefineSelectedBoxes.C
utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierEnsureCorrectRegularity.C
utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierParallelRefinement.C
utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierDistributeLeavesToProcessors.C
utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierLoadDistribution.C
utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierReduceMemoryConsumption.C
utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifierUpdateCommunicationPattern.C

utilities/octrees/meshOctree/meshOctreeInsideOutside/meshOctreeInsideOutside.C

utilities/octrees/meshOctree/meshOctreeCreator/meshOctreeCreator.C
utilities/octrees/meshOctree/meshOctreeCreator/meshOctreeCreatorAdjustOctreeToSurface.C
utilities/octrees/meshOctree/meshOctreeCreator/meshOctreeCreatorCreateOctreeBoxes.C
utilities/octrees/meshOctree/meshOctreeCreator/meshOctreeCreatorFrontalMarking.C
utilities/octrees/meshOctree/meshOctreeCreator/meshOctreeCreatorLoadDistribution.C

utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressing.C
utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressingCreation.C
utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressingGluedMesh.C
utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressingIrregularConnections.C
utilities/octrees/meshOctree/meshOctreeAddressing/meshOctreeAddressingParallelAddressing.C

utilities/octrees/meshOctree/meshOctreeAutomaticRefinement/meshOctreeAutomaticRefinement.C
utilities/octrees/meshOctree/meshOctreeAutomaticRefinement/meshOctreeAutomaticRefinementRef.C

utilities/octrees/meshOctree/meshOctree.C
utilities/octrees/meshOctree/meshOctreeCubePatches.C
utilities/octrees/meshOctree/meshOctreeNeighbourSearches.C
utilities/octrees/meshOctree/meshOctreeFindNearestSurfacePoint.C
utilities/octrees/meshOctree/meshOctreeInsideCalculations.C
utilities/octrees/meshOctree/meshOctreeParallelCommunication.C

utilities/octrees/meshOctree/refinementControls/patchRefinement/patchRefinement.C
utilities/octrees/meshOctree/refinementControls/objectRefinement/objectRefinement.C
utilities/octrees/meshOctree/refinementControls/objectRefinement/newObjectRefinement.C
utilities/octrees/meshOctree/refinementControls/objectRefinement/sphereRefinement.C
utilities/octrees/meshOctree/refinementControls/objectRefinement/lineRefinement.C
utilities/octrees/meshOctree/refinementControls/objectRefinement/coneRefinement.C
utilities/octrees/meshOctree/refinementControls/objectRefinement/boxRefinement.C
utilities/octrees/meshOctree/refinementControls/objectRefinement/hollowConeRefinement.C

utilities/triSurfaceTools/triSurface2DCheck/triSurface2DCheck.C

utilities/triSurfaceTools/triSurfaceChecks/triSurfaceChecks.C

utilities/triSurfaceTools/triSurfaceCleanupDuplicates/triSurfaceCleanupDuplicates.C
utilities/triSurfaceTools/triSurfaceCleanupDuplicates/triSurfaceCleanupDuplicatesFunctions.C

utilities/triSurfaceTools/triSurfaceCleanupDuplicateTriangles/triSurfaceCleanupDuplicateTriangles.C
utilities/triSurfaceTools/triSurfaceCleanupDuplicateTriangles/triSurfaceCleanupDuplicateTrianglesFunctions.C

utilities/triSurfaceTools/triSurfaceCopyParts/triSurfaceCopyParts.C

utilities/triSurfaceTools/triSurfacePartitioner/triSurfacePartitioner.C
utilities/triSurfaceTools/triSurfacePartitioner/triSurfacePartitionerCreateAddressing.C

utilities/triSurfaceTools/triSurfaceCurvatureEstimator/triSurfaceCurvatureEstimator.C
utilities/triSurfaceTools/triSurfaceCurvatureEstimator/triSurfaceCurvatureEstimatorCalculate.C

utilities/triSurfaceTools/triSurfaceDetectFeatureEdges/triSurfaceDetectFeatureEdges.C
utilities/triSurfaceTools/triSurfaceDetectFeatureEdges/triSurfaceDetectFeatureEdgesFunctions.C

utilities/triSurfaceTools/triSurfaceClassifyEdges/triSurfaceClassifyEdges.C
utilities/triSurfaceTools/triSurfaceClassifyEdges/triSurfaceClassifyEdgesFunctions.C

utilities/triSurfaceTools/triSurfaceImportSurfaceAsSubset/triSurfaceImportSurfaceAsSubset.C

utilities/triSurfaceTools/triSurfacePatchManipulator/triSurfacePatchManipulator.C
utilities/triSurfaceTools/triSurfacePatchManipulator/triSurfacePatchManipulatorFunctions.C

utilities/triSurfaceTools/triSurfaceRemoveFacets/triSurfaceRemoveFacets.C
utilities/triSurfaceTools/triSurfaceRemoveFacets/triSurfaceRemoveFacetsFunctions.C

utilities/triSurfaceTools/triSurfaceExtrude2DEdges/triSurfaceExtrude2DEdges.C

utilities/triSurfaceTools/triSurfaceMetaData/triSurfaceMetaData.C

cartesianMesh/cartesianMeshExtractor/cartesianMeshExtractor.C
cartesianMesh/cartesianMeshExtractor/cartesianMeshExtractorPointsAndAddressing.C
cartesianMesh/cartesianMeshExtractor/cartesianMeshExtractorPolyMesh.C
cartesianMesh/cartesianMeshExtractor/cartesianMeshExtractorDecomposeSplitHexes.C

cartesianMesh/cartesianMeshGenerator/cartesianMeshGenerator.C

cartesian2DMesh/cartesian2DMeshGenerator/cartesian2DMeshGenerator.C

tetMesh/tetMeshExtractorOctree/tetMeshExtractorOctree.C

tetMesh/tetMeshGenerator/tetMeshGenerator.C

voronoiMesh/voronoiMeshExtractor/voronoiMeshExtractor.C
voronoiMesh/voronoiMeshExtractor/voronoiMeshExtractorAddressing.C
voronoiMesh/voronoiMeshExtractor/voronoiMeshExtractorCreateMesh.C

voronoiMesh/voronoiMeshGenerator/voronoiMeshGenerator.C

utilities/dataConversion/foamToFPMA/writeMeshFPMA.C
utilities/dataConversion/foamToFPMA/fpmaMesh.C

utilities/workflowControls/workflowControls.C
)

# Building library #############################################################

# Prepend flags meant for host compiler
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if(languages MATCHES CUDA)
  set(cuda_host_flag "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>")
else()
  set(cuda_host_flag "")
endif()

# Build cfMesh library
add_library(cfmesh ${CFMESH_SRCS})
target_include_directories(cfmesh
    PUBLIC
        $<INSTALL_INTERFACE:include/cfmesh>
)
foreach(include_dir ${CFMESH_INC_DIRS})
  target_include_directories(cfmesh
      PUBLIC
          $<BUILD_INTERFACE:${include_dir}>
  )
endforeach()
target_compile_options(cfmesh
    PRIVATE
        -Wall
        -Wextra
        -Wold-style-cast
        -Wnon-virtual-dtor
        -Wno-unused-parameter
        -Wno-invalid-offsetof

    PUBLIC
        "${cuda_host_flag}-ftemplate-depth-100"
        "${cuda_host_flag}${OpenMP_CXX_FLAGS}"
)
target_compile_definitions(cfmesh
    PUBLIC
        USE_OMP
        NoSizeType
        OpenCFDSpecific
)

# Link OpenFOAM
target_link_libraries(cfmesh
    PRIVATE
        ${OPNF_LIBRARIES}
)
target_include_directories(cfmesh SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)
target_compile_definitions(cfmesh
    PUBLIC
        ${OPNF_COMPILE_DEFINITIONS}
)

# Set soname version ###########################################################

set_target_properties(cfmesh PROPERTIES VERSION ${cfmesh_VERSION}
                                        SOVERSION ${cfmesh_MAJOR_VERSION})

# Install libraries ############################################################

# Headers
install(FILES ${CFMESH_INC_FILES}
    DESTINATION include/cfmesh)

# Libraries
install(TARGETS cfmesh
    EXPORT Nemosys
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
