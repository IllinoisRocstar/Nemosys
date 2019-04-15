cmake_minimum_required(VERSION 2.8)
project(MAdLib)

message(STATUS "MAdLib Base")

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

opt(ANN "Enable ANN library" ON)

add_library(MAdLib
    Adapt/AdaptInterface.cc
    Adapt/MeshParametersManager.cc
    Adapt/constraint/Constraint.cc
    Adapt/constraint/ModelConstraintManager.cc
    Adapt/operator/DESCOp.cc
    Adapt/operator/EdgeCollapseOp.cc
    Adapt/operator/EdgeSplitOp.cc
    Adapt/operator/EdgeSwapConfig.cc
    Adapt/operator/EdgeSwapOp.cc
    Adapt/operator/FaceCollapseOp.cc
    Adapt/operator/FaceSwapOp.cc
    Adapt/operator/MAdOperatorBase.cc
    Adapt/operator/OperatorTools.cc
    Adapt/operator/RegionRemoveOp.cc
    Adapt/operator/SliverFaceHandler.cc
    Adapt/operator/SliverRegionHandler.cc
    Adapt/operator/VertexMoveOp.cc
    Adapt/output/MAdOutput.cc
    Adapt/quality/MeanRatioEvaluator.cc
    Adapt/quality/OrientedMeanRatioEvaluator.cc
    Adapt/quality/MeshQualityManager.cc
    Adapt/repositioning/GeoMatcher.cc
    Adapt/repositioning/LaplaceSmoothingOp.cc
    Adapt/repositioning/MAdElasticityOp.cc
    Adapt/repositioning/MobileObject.cc
    Adapt/sizeField/AnalyticalSField.cc
    Adapt/sizeField/AnisoMeshSize.cc
    Adapt/sizeField/BackgroundSF.cc
    Adapt/sizeField/DiscreteSF.cc
    Adapt/sizeField/IsoMeshSize.cc
    Adapt/sizeField/LocalSizeField.cc
    Adapt/sizeField/MeshSizeBase.cc
    Adapt/sizeField/NullSField.cc
    Adapt/sizeField/PWLinearSField.cc
    Adapt/sizeField/SizeFieldBase.cc
    Adapt/sizeField/SizeFieldManager.cc
    Adapt/utils/CallbackManager.cc
    Adapt/utils/History.cc
    Adapt/utils/MAdStatistics.cc
    Adapt/utils/MAdTimeManager.cc
    Adapt/utils/NodalDataManager.cc
    Adapt/utils/DistanceFunction.cc

    Common/MathUtils.cc
    Common/MAdStringFieldEvaluator.cc
    Common/MAdMessage.cc
    Common/MAdMatrix.cc
    Common/MAdMetric.cc
    Common/MAdResourceManager.cc
    Common/MAdLib.cc

    Geo/GmshModel.cc
    Geo/ModelInterface.cc
    Geo/NullModel.cc
    Geo/GmshEntities.cc
    Geo/NullEntities.cc
    Geo/MAdModel.cc

    Mesh/Balance.cc
    Mesh/CheckMesh.cc
    Mesh/CheckOrientation.cc
    Mesh/Mark.cc
    Mesh/MeshDataBase.cc
    Mesh/MeshDataBaseComm.cc
    Mesh/MeshDataBaseCommCheck.cc
    Mesh/MeshDataBaseGEntity2Physical.cc
    Mesh/MeshDataBaseInterface.cc
    Mesh/MeshDataBaseIO.cc
    Mesh/MeshDataBaseLoadBalancing.cc
    Mesh/MeshDataBaseMigration.cc
    Mesh/MeshDataBaseMiniMesh.cc
    Mesh/MeshDataBaseParallelInterface.cc
    Mesh/MeshDataBaseParallelIO.cc
    Mesh/metisAdaptiveRepart.cc
    Mesh/ParallelUtils.cc
    Mesh/PeriodicInterfaceMigration.cc
)
target_include_directories(MAdLib PRIVATE
    Adapt
    Adapt/constraint
    Adapt/operator
    Adapt/output
    Adapt/quality
    Adapt/repositioning
    Adapt/sizeField
    Adapt/utils
    Common
    Geo
    Mesh
)

target_sources(MAdLib PRIVATE Contrib/mathex/mathex.cpp)
target_include_directories(MAdLib PRIVATE Contrib/mathex)

add_subdirectory(Benchmarks/checkMesh)
add_subdirectory(Benchmarks/meshInfo)
add_subdirectory(Benchmarks/moveIt)
add_subdirectory(Benchmarks/optimize)


if(ENABLE_ANN)
  add_subdirectory(Contrib/ANN)
  target_link_libraries(MAdLib PRIVATE ANN)
  target_compile_definitions(MAdLib PRIVATE _HAVE_ANN_)
endif()
