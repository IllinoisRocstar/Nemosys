cmake_minimum_required(VERSION 3.1)

project(foamutils CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Finding libraries ############################################################

# Find OpenFOAM
# checking OpenFOAM version
string(FIND $ENV{WM_PROJECT_VERSION} "2006" OF2006)
string(FIND $ENV{WM_PROJECT_VERSION} "2012" OF2012)
if(OF2006 GREATER -1 OR OF2012 GREATER -1)
  # OpenFOAM-v2006 or OpenFOAM-v2012 present
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM
      meshTools finiteVolume dynamicMesh lagrangian
      fileFormats surfMesh distributed snappyHexMesh blockMesh
      genericPatchFields dynamicMesh topoChangerFvMesh dynamicFvMesh sampling
      fvOptions incompressibleTransportModels incompressibleTurbulenceModels
      turbulenceModels decompose decompositionMethods reconstruct meshLibrary)
else()
  message( FATAL_ERROR "Only openfoam.com versions (v2006 and v2012) 
                      are supported! CMake exiting ..." )
endif()

# Finding libraries ############################################################
find_package(OpenMP REQUIRED)

# Setting library files ########################################################

# Setting cfMesh source
set(FOAMUTILS_SRCS
  src_files/snappyMesh.C
  src_files/mergePolyMesh.C
  src_files/blockMsh.C
  src_files/createPatch.C
  src_files/foamToSurface.C
  src_files/mergeMeshes.C
  src_files/splitMeshRegions.C
  src_files/surfaceLambdaMuSmooth.C
  src_files/surfaceSplitByTopology.C
  src_files/getDicts.C
)

# Building library #############################################################

# Prepend flags meant for host compiler
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if(languages MATCHES CUDA)
  set(cuda_host_flag "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>")
else()
  set(cuda_host_flag "")
endif()

# Build FOAMUTILS library
add_library(foamutils ${FOAMUTILS_SRCS})
target_include_directories(foamutils
    PUBLIC
        $<INSTALL_INTERFACE:include/foamutils>
)
foreach(include_dir ${CMAKE_CURRENT_LIST_DIR}/include)
  target_include_directories(foamutils
      PUBLIC
          $<BUILD_INTERFACE:${include_dir}>
  )
endforeach()
target_compile_options(foamutils
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
target_compile_definitions(foamutils
    PUBLIC
        USE_OMP
        NoSizeType
        OpenCFDSpecific
)

# Link OpenFOAM
target_link_libraries(foamutils
    PRIVATE
        ${OPNF_LIBRARIES}
)
target_include_directories(foamutils SYSTEM
    PRIVATE
        ${OPNF_INCLUDE_DIRS}
)

target_compile_definitions(foamutils
    PUBLIC
        ${OPNF_COMPILE_DEFINITIONS}
)

# Install libraries ############################################################

# Headers
install(FILES ${CFMESH_INC_FILES}
    DESTINATION ${CMAKE_CURRENT_LIST_DIR}/include)

# Libraries
install(TARGETS foamutils
    EXPORT Nemosys
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
