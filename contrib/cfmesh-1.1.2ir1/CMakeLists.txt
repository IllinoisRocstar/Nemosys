cmake_minimum_required(VERSION 3.1)

project(cfmesh CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set project metadata #########################################################

set(cfmesh_MAJOR_VERSION 1)
set(cfmesh_MINOR_VERSION 1)
set(cfmesh_PATCH_VERSION 2)
set(cfmesh_EXTRA_VERSION "ir1")
set(cfmesh_VERSION "${cfmesh_MAJOR_VERSION}.${cfmesh_MINOR_VERSION}")
set(cfmesh_VERSION
    "${cfmesh_VERSION}.${cfmesh_PATCH_VERSION}${cfmesh_EXTRA_VERSION}")

# Finding libraries ############################################################

# Find OpenFOAM
# checking OpenFOAM version
string(FIND $ENV{WM_PROJECT_VERSION} "4" OF4)
if(OF4 GREATER -1)
  # OpenFOAM 4.x
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM triSurface
      meshTools finiteVolume foamToVTK dynamicMesh lagrangian
      fileFormats surfMesh triSurface distributed snappyHexMesh blockMesh
      genericPatchFields
      edgeMesh)
else()
  # 5.x and 6.x series
  find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM triSurface
      meshTools finiteVolume foamToVTK dynamicMesh lagrangian
      fileFormats surfMesh triSurface distributed snappyHexMesh blockMesh
      genericPatchFields)
endif()

# Add sub-projects #############################################################

add_subdirectory(meshLibrary)
#add_subdirectory(executables)
#add_subdirectory(utilities)
