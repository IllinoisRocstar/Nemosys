cmake_minimum_required(VERSION 3.1)
project(MAdLib CXX)
set(CMAKE_CXX_STANDARD 11)

set(MAdLib_MAJOR_VERSION 1)
set(MAdLib_MINOR_VERSION 3)
set(MAdLib_PATCH_VERSION 0)
set(MAdLib_EXTRA_VERSION "")
set(MAdLib_VERSION "${MAdLib_MAJOR_VERSION}.${MAdLib_MINOR_VERSION}")
set(MAdLib_VERSION
    "${MAdLib_VERSION}.${MAdLib_PATCH_VERSION}${MAdLib_EXTRA_VERSION}")

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

opt(ANN "Enable ANN library" ON)
opt(GMSH "Enable Gmsh geometric model" ON)
opt(BENCHMARKS "Enable Benchmark executables" OFF)

# Check for <unistd.h> header
include(CheckIncludeFile)
check_include_file(unistd.h HAVE_UNISTD_H)

# Set MAdLib include files
file(GLOB_RECURSE Adapt_FILES Adapt/*.h)
file(GLOB Common_FILES Common/*.h)
file(GLOB Geo_FILES Geo/*.h)
file(GLOB Mesh_FILES Mesh/*.h)
set(MADLIB_INC_FILES ${Adapt_FILES} ${Common_FILES} ${Geo_FILES} ${Mesh_FILES})

# Set MAdLib source files
set(MADLIB_SRCS
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

# Build MAdLib library
add_library(MAdLib ${MADLIB_SRCS})

# Set MAdLib include directory
target_include_directories(MAdLib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt/constraint>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt/operator>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt/output>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt/quality>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt/repositioning>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt/sizeField>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Adapt/utils>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Common>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geo>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Mesh>

        $<INSTALL_INTERFACE:include/MAdLib>
)

# Set MAdLib compiler definitions
if(HAVE_UNISTD_H)
  target_compile_definitions(MAdLib PUBLIC HAVE_UNISTD_H)
endif()

# Link mathex library
target_sources(MAdLib PRIVATE Contrib/mathex/mathex.cpp)
target_include_directories(MAdLib PRIVATE Contrib/mathex)

# Generate export header #######################################################
include(GenerateExportHeader)
generate_export_header(MAdLib)
list(APPEND ${Nemosys_INC_FILES} madlib_export.h)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/madlib_export.h
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/Common)

# Link ANN library
if(ENABLE_ANN)
  add_subdirectory(Contrib/ANN)
  target_link_libraries(MAdLib PUBLIC ANN)
  target_compile_definitions(MAdLib PUBLIC _HAVE_ANN_)
endif()

# Link Gmsh library
if(ENABLE_GMSH)
  find_library(GMSH_LIB gmsh HINTS lib)
  find_file(GMSH_HDR gmsh.h HINTS include)
  get_filename_component(GMSH_INCPATH ${GMSH_HDR} PATH)
  message(STATUS "gmsh include path ${GMSH_INCPATH}")
  message(STATUS "gmsh libs ${GMSH_LIB}")

#  find_package(OpenCASCADE REQUIRED) # Used by Gmsh

  target_link_libraries(MAdLib PRIVATE ${GMSH_LIB})
  target_include_directories(MAdLib SYSTEM PRIVATE ${GMSH_INCPATH})
  target_compile_definitions(MAdLib PUBLIC _HAVE_GMSH_)
endif()

# Add benchmark executables
if(ENABLE_BENCHMARKS)
  add_subdirectory(Benchmarks/checkMesh)
  add_subdirectory(Benchmarks/meshInfo)
  add_subdirectory(Benchmarks/moveIt)
  add_subdirectory(Benchmarks/optimize)
endif()

# Set soname version
set_target_properties(MAdLib PROPERTIES VERSION ${MAdLib_VERSION}
                                        SOVERSION ${MAdLib_MAJOR_VERSION})

# Install MAdLib headers
install(FILES ${MADLIB_INC_FILES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MAdLib)

# Install MAdLib library
install(TARGETS MAdLib
    EXPORT Nemosys
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
