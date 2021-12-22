# TODO: Make this hierarchical
cmake_minimum_required(VERSION 3.9)

# Set project metadata #########################################################

set(NEMOSYS_MAJOR_VERSION 0)
set(NEMOSYS_MINOR_VERSION 62)
set(NEMOSYS_PATCH_VERSION 1)
set(NEMOSYS_EXTRA_VERSION "" CACHE STRING "NEMoSys extra version string")
set(NEMOSYS_VERSION "${NEMOSYS_MAJOR_VERSION}.${NEMOSYS_MINOR_VERSION}")
set(NEMOSYS_VERSION
    "${NEMOSYS_VERSION}.${NEMOSYS_PATCH_VERSION}${NEMOSYS_EXTRA_VERSION}")
set(NEMOSYS_SHORT_LICENSE "GNU General Public License v2 or later (GPLv2+)")
set(NEMOSYS_HOMEPAGE_URL "http://promesh.com")
set(NEMOSYS_AUTHOR
    "Alessandro Gondolo"
    "Andrew Fu"
    "Sachin Natesh"
    "Akash Patel"
    "Seth Pemberton"
    "Masoud Safdari"
    "Pavle Simonovic"
    "Taiyo Wilson")
set(NEMOSYS_AUTHOR_EMAIL "promesh@illinoisrocstar.com")
set(NEMOSYS_PYTHON_CLASSIFIERS
    "Development Status :: 2 - Pre-Alpha"
    "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)"
    "Natural Language :: English"
    "Operating System :: Microsoft :: Windows :: Windows 10"
    "Operating System :: POSIX"
    "Programming Language :: C"
    "Programming Language :: C++"
    "Topic :: Scientific/Engineering")
set(NEMOSYS_DESCRIPTION "The Nuclear Energy Modeling System is a modular, \
extensible resource designed for use in typical application development \
systems as well as distributed web-services environments. The project focus is \
providing a framework for robust, automated mesh generation, mesh quality \
analysis, adaptive mesh refinement, and data transfer between arbitrary \
meshes.")
set(PROJECT_NAME "NEMoSys version ${NEMOSYS_VERSION}")
project(${PROJECT_NAME} CXX C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set project defaults #########################################################

# Set a default build type if none was specified
set(default_build_type "Release")
# TODO: Debug build fails on Windows. Figure out why and fix.
if(EXISTS "${CMAKE_SOURCE_DIR}/.git" AND NOT WIN32)
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
      "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Setting RPATH. Must be done before any add_<target>() commands.
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Installation Directories
include(GNUInstallDirs)

# Setting options ##############################################################

# package module
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

# Macros and functions
include(cmake/NemosysUtils.cmake)

# This variable controls the default value of the options, which are normally
# set to ON (useful if you want to configure a minimal version of NEMoSys: e.g.,
#  "cmake -DDEFAULT=0 -DENABLE_BUILD_UTILS=1")
set(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")
# creating code compilation database
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_BUILD_TYPE_OPTIONS
    "None" "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

# set to ${DEFAULT} to use by DEFAULT
option(ENABLE_MPI "Enable MPI support" OFF)
option(ENABLE_TESTING "Enable testing" ${DEFAULT})
option(ENABLE_BUILD_UTILS "Build utilities" OFF)
option(ENABLE_PYTHON_BINDINGS "Enable Python bindings" OFF)
option(ENABLE_CONSRV_SURFACE_TRANSFER "Enable Conservative Surface Transfer" OFF)
option(ENABLE_CONSRV_VOLUME_TRANSFER "Enable Conservative Volume Transfer" OFF)
option(ENABLE_CGNS "Enable CGNS extensions (Used for Rocstar utilities)" OFF)
option(ENABLE_EPIC "Enable EPIC preprocessor" OFF)
option(ENABLE_HDF5 "ENABLE HDF5 extensions (Used by some CGNS and Netgen builds)" OFF)
option(ENABLE_METIS "Enable METIS partitioner" ${DEFAULT})
option(ENABLE_NETGEN "Enable Netgen meshing engine" ${DEFAULT})
option(ENABLE_OMEGAH_CUDA "Enable GPU acceleration for Omega_h (Requires Kokkos built \
with CUDA and linked dynamically)" OFF)
option(ENABLE_OPENCASCADE "Enable OpenCASCADE kernel support (Used by some GMSH builds)" ${DEFAULT})
option(ENABLE_MLAMR "Enable machine learning module in AMR" OFF)
option(ENABLE_OPENMP "Enable OpenMP" ${DEFAULT})

option(ENABLE_GPL "Allow linking against GPL libraries, making combined project fall under GPL license requirements" OFF)
option(ENABLE_CFMSH "Enable cfMesh Meshing engine" OFF)
option(ENABLE_GMSH "Enable GMSH" OFF)
option(ENABLE_TEMPLATE_MESH "Enable meshing templates (requires gmsh)" OFF)


# Check for option interdependencies ###########################################

if(ENABLE_MLAMR)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(WIN32)
  if(ENABLE_CONSRV_SURFACE_TRANSFER)
    message(FATAL_ERROR "Conservative surface transfer is not supported on Windows.")
  endif()
  if(ENABLE_CONSRV_VOLUME_TRANSFER)
    message(FATAL_ERROR "Conservative volume transfer is not supported on Windows.")
  endif()
endif()

if(ENABLE_CONSRV_VOLUME_TRANSFER AND NOT ENABLE_MPI)
  message(FATAL_ERROR "Conservative volume transfer requires MPI.")
endif()

if(ENABLE_OMEGAH_CUDA)
  if(MSVC)
    message(FATAL_ERROR "MSVC and CUDA not supported")
  endif()
endif()

if(ENABLE_CFMSH AND NOT ENABLE_OPENMP)
  message(FATAL_ERROR "cfMesh requires OpenMP.")
endif()

if(NOT ENABLE_GPL)
  if(ENABLE_GMSH)
    message(FATAL_ERROR "Linking against Gmsh requires ENABLE_GPL flag to be set to ON")
  endif()
  if(ENABLE_CFMSH)
    message(FATAL_ERROR "Linking against cfMesh requires ENABLE_GPL flag to be set to ON")
  endif()
endif()

if(ENABLE_TEMPLATE_MESH AND NOT ENABLE_GMSH)
  message(FATAL_ERROR "Template mesh requires ENABLE_GMSH to be set to ON")
endif()

# Check system configuration ###################################################

# Check for system-depended headers
include(CheckIncludeFile)
check_include_file(glob.h HAVE_GLOB_H)

# Finding libraries ############################################################

# Find MPI
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  # include mpi headers
  include_directories(${MPI_INCLUDE_PATH})
  if(MPI_LIBRARY)
    set(MPI_LIBRARIES ${MPI_LIBRARY})
  endif(MPI_LIBRARY)
  if(MPI_EXTRA_LIBRARY)
    set(MPI_LIBRARIES ${MPI_LIBRARIES} "${MPI_EXTRA_LIBRARY}")
  endif(MPI_EXTRA_LIBRARY)
  # Needed for mpich 2
  add_definitions("-DMPICH_IGNORE_CXX_SEEK")
  set(Omega_h_USE_MPI OFF)
endif()

# Find VTK
# v8.2.0 required for thread-safe point and cell locators
find_package(VTK 8.2.0 REQUIRED)

# UseVTK.cmake file does not correctly set VTK headers as external.
#include(${VTK_USE_FILE})
set(VTK_LIB_DIR ${vtksys_RUNTIME_LIBRARY_DIRS})
list(APPEND LIBRARY_DIRS ${VTK_LIB_DIR})
message(STATUS "vtk include path: ${VTK_INCLUDE_DIRS}")
message(STATUS "vtk libs ${VTK_LIB_DIR}")

# Find ExodusII
find_library(EXO_LIB NAMES exoIIv2 exodus exoIIv2c)
find_file(EXO_HDR exodusII.h)
get_filename_component(EXO_LIB_PATH ${EXO_LIB} PATH)
get_filename_component(EXO_INCPATH ${EXO_HDR} PATH)
list(APPEND LIBRARY_DIRS ${EXO_LIB_PATH})
message(STATUS "exodus include path ${EXO_INCPATH}")
message(STATUS "exodus libs ${EXO_LIB}")

# Find NetCDF
if(WIN32)
  find_package(netCDF REQUIRED)
  list(APPEND LIBRARY_DIRS ${netCDF_LIB_DIR} ${netCDF_LIB_DIR}/../bin)
  message(STATUS "netcdf include path ${netCDF_INCLUDE_DIR}")
  message(STATUS "netcdf libs ${netCDF_LIB_DIR}")
endif()

# Find OpenCascade
if(ENABLE_OPENCASCADE)
  get_directory_property(PRE_OCC_COMPILE_DEFS COMPILE_DEFINITIONS)
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "")
  find_package(OpenCASCADE REQUIRED)
  list(APPEND LIBRARY_DIRS ${OpenCASCADE_LIBRARY_DIR})
  message(STATUS "OpenCASCADE include path ${OpenCASCADE_INCLUDE_DIR}")
  message(STATUS "OpenCASCADE libs ${OpenCASCADE_LIBRARY_DIR}")
  get_directory_property(OCC_COMPILE_DEFINITIONS COMPILE_DEFINITIONS)
  set_directory_properties(PROPERTIES COMPILE_DEFINITIONS "${PRE_OCC_COMPILE_DEFS}")
endif()

# Find Gmsh
if(ENABLE_GMSH)
  find_package(Gmsh REQUIRED)
  get_filename_component(GMSH_LIB_PATH ${Gmsh_LIBRARIES} PATH)
  list(APPEND LIBRARY_DIRS ${GMSH_LIB_PATH})
  message(STATUS "gmsh include path ${Gmsh_INCLUDE_DIRS}")
  message(STATUS "gmsh libs ${Gmsh_LIBRARIES}")
endif()

# Find CGNS
if(ENABLE_CGNS)
  if(WIN32 AND BUILD_SHARED_LIBS)
    find_library(CGNS_LIB cgnsdll)
  else()
    find_library(CGNS_LIB cgns)
  endif()
  find_file(CGNS_HDR cgnslib.h)
  get_filename_component(CGNS_INCPATH ${CGNS_HDR} PATH)
  get_filename_component(CGNS_LIB_PATH ${CGNS_LIB} PATH)
  list(APPEND LIBRARY_DIRS ${CGNS_LIB_PATH} ${CGNS_LIB_PATH}/../bin)
  message(STATUS "cgns include path ${CGNS_INCPATH}")
  message(STATUS "cgns libs ${CGNS_LIB}")
  set(HAVE_CGNS TRUE)
endif()

# Find HDF5
if(ENABLE_HDF5)
  if(WIN32) # Test new HDF5 logic on UNIX
    find_package(hdf5 REQUIRED COMPONENTS CXX CXX_HL CONFIG)
    NEM_get_imported_location(HDF5_LIBRARIES hdf5::hdf5-shared HDF5_LIB_PATH)
    list(APPEND LIBRARY_DIRS ${HDF5_LIB_PATH})
    message(STATUS "hdf5 include path ${HDF5_INCLUDE_DIR}")
    message(STATUS "hdf5 libs ${HDF5_LIBRARIES}")
    set(HAVE_HDF5 TRUE)
    set(HAVE_HDF5_CONFIG TRUE)
  else()
    find_package(HDF5 REQUIRED COMPONENTS CXX)
    list(GET ${HDF5_LIBRARIES} 0 HDF5_LIBRARY)
    get_filename_component(HDF5_LIB_PATH ${HDF5_LIBRARY} PATH)
    list(APPEND LIBRARY_DIRS ${HDF5_LIB_PATH})
    message(STATUS "hdf5 include path ${HDF5_INCLUDE_DIRS}")
    message(STATUS "hdf5 libs ${HDF5_LIBRARIES}")
    set(HAVE_HDF5 TRUE)
    set(HAVE_HDF5_CONFIG FALSE)
  endif()
endif()

# Find METIS
if(ENABLE_METIS)
  find_library(METIS_LIB metis)
  find_file(METIS_HDR metis.h)
  get_filename_component(METIS_INCPATH ${METIS_HDR} PATH)
  get_filename_component(METIS_LIB_PATH ${METIS_LIB} PATH)
  list(APPEND LIBRARY_DIRS ${METIS_LIB_PATH})
  message(STATUS "metis include path ${METIS_INCPATH}")
  message(STATUS "metis libs ${METIS_LIB}")
  set(HAVE_METIS TRUE)
endif()

# Find Netgen
if(ENABLE_NETGEN)
  find_package(Netgen)
  if(NOT Netgen_FOUND)
    find_library(NETGEN_LIB nglib)
    get_filename_component(NETGEN_LIBRARY_DIR ${NETGEN_LIB} PATH)
    foreach(comp mesh visual csg geom2d occ occvis stl interface)
      find_library(NETGEN_LIB_${comp} ${comp} HINTS ${NETGEN_LIBRARY_DIR})
      list(APPEND NETGEN_LIB ${NETGEN_LIB_${comp}})
    endforeach()

    find_file(NETGEN_HDR nglib.h PATH /usr/include/netgen-mesher)
    get_filename_component(NETGEN_INCLUDE_DIRS ${NETGEN_HDR} PATH)

    set(HAVE_NETGEN_OLD TRUE)
  endif()
  list(APPEND LIBRARY_DIRS ${NETGEN_LIBRARY_DIR} ${NETGEN_LIBRARY_DIR}/../bin)
  message(STATUS "netgen include path ${NETGEN_INCLUDE_DIRS}")
  message(STATUS "netgen libs ${NETGEN_LIBRARY_DIR}")

  # Find zlib
  if(WIN32)
    add_definitions(-D_AFXDLL) # Needed by Netgen headers.
    find_package(ZLIB REQUIRED) # Find zlib headers since Netgen doesn't
    include_directories(${ZLIB_INCLUDE_DIRS})
    get_filename_component(ZLIB_LIB_PATH ${ZLIB_LIBRARIES} PATH)
    list(APPEND LIBRARY_DIRS ${ZLIB_LIB_PATH} ${ZLIB_LIB_PATH}/../bin)
    message(STATUS "zlib include path ${ZLIB_INCLUDE_DIRS}")
    message(STATUS "zlib libs ${ZLIB_LIBRARIES}")
  endif()
endif()

# Find boost
# Needed if there is no glob.h on the system
if(ENABLE_CFMSH OR (ENABLE_BUILD_UTILS AND HAVE_CGNS AND NOT HAVE_GLOB_H))
  if(WIN32)
    set(Boost_USE_STATIC_LIBS ON) # needed ON on WIN32.
  endif()
  set(Boost_USE_MULTITHREADED ON)
  find_package(Boost 1.45.0 REQUIRED COMPONENTS filesystem)
  if(ENABLE_CFMSH)
    list(APPEND LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
  endif()
endif()

# Find cfMesh and OpenFOAM
if(ENABLE_CFMSH)
  # checking OpenFOAM version
  string(FIND $ENV{WM_PROJECT_VERSION} "2006" OF2006)
  if(OF2006 GREATER -1)
    # OpenFOAM-v2006 present
    find_package(OpenFOAM REQUIRED COMPONENTS OpenFOAM
        meshTools finiteVolume dynamicMesh lagrangian
        fileFormats surfMesh distributed snappyHexMesh blockMesh
        genericPatchFields dynamicMesh topoChangerFvMesh dynamicFvMesh sampling
        fvOptions incompressibleTransportModels incompressibleTurbulenceModels
        turbulenceModels decompose decompositionMethods reconstruct meshLibrary)
  else()
    message(FATAL_ERROR "Only openfoam.com v2006 is supported.")
  endif()

  message(STATUS "OpenFOAM include directory ${OPNF_INCLUDE_DIRS}")
  message(STATUS "OpenFOAM link directory ${OPNF_LIBRARY_DIRS}")

  list(APPEND LIBRARY_DIRS ${OPNF_LIBRARY_DIRS})
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

if(ENABLE_CONSRV_SURFACE_TRANSFER)
  find_package(IMPACT REQUIRED)
endif()

# Setting compile and link options #############################################

# Adding flags to Debug
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# Compiled objects output paths
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Find contrib libraries #######################################################

# Note check for enabling CUDA before any other contrib libraries are added
# so that they can check ENABLED_LANGUAGES and use generator expressions to be
# sure C++ flags are passed to the host compiler and not nvcc.
if(ENABLE_OMEGAH_CUDA)
  set(Omega_h_USE_CUDA ON CACHE INTERNAL "")
  set(Omega_h_USE_Kokkos ON CACHE INTERNAL "")
  enable_language(CUDA)
else()
  set(Omega_h_USE_CUDA OFF CACHE INTERNAL "")
  set(Omega_h_USE_Kokkos OFF CACHE INTERNAL "")
endif()
add_subdirectory(contrib/omega_h-9.31.2)
list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/omega_h-9.31.2/src)

# Find ANN
add_subdirectory(contrib/ann-1.1.2)
list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/ann-1.1.2)

# Find MAdLib (make sure done after ANN)
if(ENABLE_GMSH)
  add_subdirectory(contrib/madlib-1.3.0) # includes ANN
  list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0)
  list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0/Contrib/ANN)
endif()

# Find jsoncons
add_subdirectory(contrib/jsoncons-0.159.0)

# Find Eigen
add_subdirectory(contrib/eigen-3.3.7)

if(ENABLE_CONSRV_VOLUME_TRANSFER)
  add_subdirectory(contrib/libsupermesh-1.0.1)
  # TODO: Add LIBRARY_DIRS
endif()

if(ENABLE_EPIC)
  add_subdirectory(contrib/quickhull-4f65e08)
  add_subdirectory(contrib/kmeans-5b9691d)
  list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/kmeans-5b9691d)
  add_subdirectory(contrib/fast-cpp-csv-parser-75600d0)
endif()

if(ENABLE_OPENCASCADE)
  add_subdirectory(contrib/smesh-9.6.0)
  list(APPEND LIBRARY_DIRS
          ${CMAKE_BINARY_DIR}/contrib/smesh-9.6.0/Kernel
          ${CMAKE_BINARY_DIR}/contrib/smesh-9.6.0/Geom
          ${CMAKE_BINARY_DIR}/contrib/smesh-9.6.0/SMESH)
  if(ENABLE_NETGEN)
    list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/smesh-9.6.0/NETGENPlugin)
  endif()
endif()

if(ENABLE_CFMSH)
  # Find foamutils
  add_subdirectory(contrib/foamutils)

  # Setting include files for CfMesh
  file(GLOB cfmsh_H_files contrib/cfmesh_headers/*.H)
  file(GLOB cfmsh_C_files contrib/cfmesh_headers/*.C)
  set(CFMSH_INC_FILES ${cfmsh_H_files} ${cfmsh_C_files})
endif()

# Setting library files ########################################################

# Setting Nemosys source
# List both headers and source files. Listed headers that are in include/ will
# be installed. Private headers should be in src/
set(NEMOSYS_SRCS
        src/Drivers/NemJsonMacros.H
        src/Drivers/DriverJsonTypeTraits.H
        include/Drivers/NemDriver.H
        src/Drivers/NemDriver.C

        src/Drivers/AutoVerificationJson.H
        include/Drivers/AutoVerificationDriver.H
        src/Drivers/AutoVerificationDriver.C

        src/Drivers/Conversion/ConversionJson.H
        include/Drivers/Conversion/ConversionDriver.H
        src/Drivers/Conversion/ConversionDriver.C
        include/Drivers/Conversion/GmshToExoConversionDriver.H
        src/Drivers/Conversion/GmshToExoConversionDriver.C
        include/Drivers/Conversion/GmshToVtkConversionDriver.H
        src/Drivers/Conversion/GmshToVtkConversionDriver.C
        include/Drivers/Conversion/ManipExoConversionDriver.H
        src/Drivers/Conversion/ManipExoConversionDriver.C
        include/Drivers/Conversion/SmartConversionDriver.H
        src/Drivers/Conversion/SmartConversionDriver.C
        include/Drivers/Conversion/VtkHexToTetConversionDriver.H
        src/Drivers/Conversion/VtkHexToTetConversionDriver.C
        include/Drivers/Conversion/VtkToCobaltConversionDriver.H
        src/Drivers/Conversion/VtkToCobaltConversionDriver.C
        include/Drivers/Conversion/VtkToPatranConversionDriver.H
        src/Drivers/Conversion/VtkToPatranConversionDriver.C
        include/Drivers/Conversion/VtkToPntConversionDriver.H
        src/Drivers/Conversion/VtkToPntConversionDriver.C

        src/Drivers/MeshGen/MeshGenJson.H
        include/Drivers/MeshGen/MeshGenDriver.H
        src/Drivers/MeshGen/MeshGenDriver.C
        src/Drivers/MeshQualityJson.H
        include/Drivers/MeshQualityDriver.H
        src/Drivers/MeshQualityDriver.C
        src/Drivers/Refine/RefineJson.H
        include/Drivers/Refine/RefineDriver.H
        src/Drivers/Refine/RefineDriver.C
        include/Drivers/Refine/OmegahRefineDriver.H
        src/Drivers/Refine/OmegahRefineDriver.C
        include/Drivers/TransferDriver.H
        src/Drivers/TransferJson.H
        src/Drivers/TransferDriver.C

        src/Drivers/PackMesh/PackMeshJson.H
        include/Drivers/PackMesh/PackMeshDriver.H
        src/Drivers/PackMesh/PackMeshDriver.C

        include/Geometry/rocPackShape.H
        src/Geometry/rocPackShape.C
        include/Geometry/hmxShape.H
        src/Geometry/hmxShape.C
        include/Geometry/petnShape.H
        src/Geometry/petnShape.C
        include/Geometry/icosidodecahedronShape.H
        src/Geometry/icosidodecahedronShape.C

        include/InputGeneration/inputGen.H
        src/InputGeneration/inputGen.C

        include/Integration/Cubature.H
        src/Integration/Cubature.C

        include/Mesh/meshBase.H
        src/Mesh/meshBase.C
        include/Mesh/geoMeshBase.H
        src/Mesh/geoMeshBase.C
        src/Mesh/geoMeshBaseReconstructGeo.C
        include/Mesh/geoMeshFactory.H
        src/Mesh/geoMeshFactory.C
        include/Mesh/cobalt.H
        src/Mesh/cobalt.C
        include/Mesh/diffMesh.H
        src/Mesh/diffMesh.C
        include/Mesh/exoMesh.H
        src/Mesh/exoMesh.C
        include/Mesh/exoGeoMesh.H
        src/Mesh/exoGeoMesh.C
        include/Mesh/inpGeoMesh.H
        src/Mesh/inpGeoMesh.C
        include/Mesh/oshGeoMesh.H
        src/Mesh/oshGeoMesh.C
        include/Mesh/patran.H
        src/Mesh/patran.C
        include/Mesh/pntMesh.H
        src/Mesh/pntMesh.C
        include/Mesh/vtkGeoMesh.H
        src/Mesh/vtkGeoMesh.C
        include/Mesh/vtkMesh.H
        src/Mesh/vtkMesh.C

        include/MeshGeneration/meshingParams.H
        include/MeshGeneration/meshGen.H

        src/MeshOperation/dataSetRegionBoundaryFilter.H
        src/MeshOperation/dataSetRegionBoundaryFilter.C
        src/MeshOperation/mergeCells.H
        src/MeshOperation/mergeCells.C
        include/MeshOperation/meshSrch.H
        src/MeshOperation/meshSrch.C

        include/MeshPartitioning/meshPartitioner.H
        src/MeshPartitioning/meshPartitioner.C

        include/MeshQuality/MeshQuality.H
        src/MeshQuality/MeshQuality.C
        include/MeshQuality/meshQualityParams.H

        include/PatchRecovery/orthoPoly1D.H
        src/PatchRecovery/orthoPoly1D.C
        include/PatchRecovery/orthoPoly3D.H
        src/PatchRecovery/orthoPoly3D.C
        include/PatchRecovery/patchRecovery.H
        src/PatchRecovery/patchRecovery.C
        include/PatchRecovery/polyApprox.H
        src/PatchRecovery/polyApprox.C

        include/Services/srvBase.H
        src/Services/srvBase.C
        include/Services/omegahRefineSrv.H
        src/Services/omegahRefineSrv.C

        include/SizeFieldGeneration/SizeFieldGen.H
        src/SizeFieldGeneration/SizeFieldGen.C
        include/SizeFieldGeneration/GradSizeField.H
        src/SizeFieldGeneration/GradSizeField.C
        include/SizeFieldGeneration/ValSizeField.H
        src/SizeFieldGeneration/ValSizeField.C
        include/SizeFieldGeneration/Z2ErrorSizeField.H
        src/SizeFieldGeneration/Z2ErrorSizeField.C

        include/SolutionVerification/OrderOfAccuracy.H
        src/SolutionVerification/OrderOfAccuracy.C

        include/Transfer/TransferBase.H
        src/Transfer/TransferBase.C
        include/Transfer/FETransfer.H
        src/Transfer/FETransfer.C

        src/AuxiliaryFunctions.H

        )

# Optional source files
if(HAVE_CGNS)
  list(APPEND NEMOSYS_SRCS
          src/Drivers/RocPartCommGenJson.H
          include/Drivers/RocPartCommGenDriver.H
          src/Drivers/RocPartCommGenDriver.C

          include/IO/cgnsWriter.H
          src/IO/cgnsWriter.C
          )
  if(ENABLE_GMSH)
    list(APPEND NEMOSYS_SRCS
            include/IO/cgnsAnalyzer.H
            src/IO/cgnsAnalyzer.C
            include/IO/rocstarCgns.H
            src/IO/rocstarCgns.C

            include/MeshPartitioning/meshStitcher.H
            src/MeshPartitioning/meshStitcher.C
            )
  endif()
endif()


if(ENABLE_HDF5)
  list(APPEND NEMOSYS_SRCS
          include/IO/hdf5Reader.H
          src/IO/hdf5Reader.C
          )
  if(ENABLE_GMSH)
    list(APPEND NEMOSYS_SRCS
            src/Drivers/ProteusJson.H
            include/Drivers/ProteusDriver.H
            src/Drivers/ProteusDriver.C
            include/IO/proteusHdf5.H
            src/IO/proteusHdf5.C
            )
  endif()
endif()


if(ENABLE_EPIC)
  list(APPEND NEMOSYS_SRCS
      src/Drivers/InputGenJson.H
      include/Drivers/InputGenDriver.H
      src/Drivers/InputGenDriver.C
      include/Geometry/convexContainer.H
      src/Geometry/convexContainer.C
      include/InputGeneration/ep16Post.H
      src/InputGeneration/ep16Post.C
      include/InputGeneration/ep16Prep.H
      src/InputGeneration/ep16Prep.C
  )
endif()

if(ENABLE_NETGEN)
  list(APPEND NEMOSYS_SRCS
          include/Drivers/MeshGen/NetgenMeshGenDriver.H
          src/Drivers/MeshGen/NetgenMeshGenDriver.C

          include/MeshGeneration/netgenParams.H
          include/MeshGeneration/netgenGen.H
          src/MeshGeneration/netgenGen.C
          )
endif()

if(ENABLE_TEMPLATE_MESH)
  list(APPEND NEMOSYS_SRCS
      include/Drivers/TemplateMeshDriver.H
      src/Drivers/TemplateMeshJson.H
      src/Drivers/TemplateMeshDriver.C
  )
endif()

if(ENABLE_CFMSH)
  list(APPEND NEMOSYS_SRCS
      include/Drivers/Conversion/FoamToVtkConversionDriver.H
      src/Drivers/Conversion/FoamToVtkConversionDriver.C
      include/Drivers/Conversion/VtkToFoamConversionDriver.H
      src/Drivers/Conversion/VtkToFoamConversionDriver.C
      include/Drivers/MeshGen/BlockMeshMeshGenDriver.H
      src/Drivers/MeshGen/BlockMeshMeshGenDriver.C
      include/Drivers/MeshGen/CFMeshMeshGenDriver.H
      src/Drivers/MeshGen/CFMeshMeshGenDriver.C
      include/Drivers/MeshGen/SnappyMeshMeshGenDriver.H
      src/Drivers/MeshGen/SnappyMeshMeshGenDriver.C
      include/Drivers/Refine/FoamRefineDriver.H
      src/Drivers/Refine/FoamRefineDriver.C
      include/Drivers/PackMesh/HexPackMeshDriver.H
      src/Drivers/PackMesh/HexPackMeshDriver.C

      include/Mesh/foamMesh.H
      src/Mesh/foamMesh.C
      include/Mesh/foamGeoMesh.H
      src/Mesh/foamGeoMesh.C

      include/MeshGeneration/blockMeshGen.H
      src/MeshGeneration/blockMeshGen.C
      include/MeshGeneration/blockMeshParams.H
      include/MeshGeneration/cfmeshGen.H
      src/MeshGeneration/cfmeshGen.C
      include/MeshGeneration/cfmeshParams.H
      include/MeshGeneration/snappymeshGen.H
      src/MeshGeneration/snappymeshGen.C
      include/MeshGeneration/snappymeshParams.H

      include/MeshManipulationFoam/MeshManipulationFoam.H
      src/MeshManipulationFoam/MeshManipulationFoam.C
      include/MeshManipulationFoam/MeshManipulationFoamParams.H

      include/MeshQuality/cfmeshQualityParams.H

      include/Refinement/AMRFoam.H
      src/Refinement/AMRFoam.C
  )
  if(ENABLE_GMSH)
    list(APPEND NEMOSYS_SRCS
            include/Drivers/Conversion/FoamToMshConversionDriver.H
            src/Drivers/Conversion/FoamToMshConversionDriver.C
            )
  endif()
endif()

if(ENABLE_CONSRV_VOLUME_TRANSFER)
  list(APPEND NEMOSYS_SRCS
      include/Transfer/ConservativeVolumeTransfer.H
      src/Transfer/ConservativeVolumeTransfer.C
  )
endif()

if(ENABLE_CONSRV_SURFACE_TRANSFER)
  list(APPEND NEMOSYS_SRCS
      include/Transfer/ConservativeSurfaceTransfer.H
      src/Transfer/ConservativeSurfaceTransfer.C
  )
endif()

if(ENABLE_GMSH)
  list(APPEND NEMOSYS_SRCS
      include/Drivers/MeshGen/GmshMeshGenDriver.H
      src/Drivers/MeshGen/GmshMeshGenDriver.C
      include/Drivers/PackMesh/SurfacePackMeshDriver.H
      src/Drivers/PackMesh/SurfacePackMeshDriver.C
      include/Drivers/Refine/SizeFieldRefineDriver.H
      src/Drivers/Refine/SizeFieldRefineDriver.C
      include/Drivers/Refine/UniformRefineDriver.H
      src/Drivers/Refine/UniformRefineDriver.C
      include/Drivers/Refine/Z2RefineDriver.H
      src/Drivers/Refine/Z2RefineDriver.C

      include/Geometry/rocPack.H
      src/Geometry/rocPack.C

      src/Mesh/gmshTypes.H
      include/Mesh/gmshMesh.H
      src/Mesh/gmshMesh.C
      include/Mesh/gmshGeoMesh.H
      src/Mesh/gmshGeoMesh.C

      include/MeshGeneration/gmshGen.H
      include/MeshGeneration/gmshParams.H
      src/MeshGeneration/gmshGen.C

      include/Refinement/Refine.H
      src/Refinement/Refine.C
      )
endif()


if(ENABLE_OPENCASCADE)
  list(APPEND NEMOSYS_SRCS
      src/Drivers/NucMeshJson.H
      include/Drivers/NucMeshDriver.H
      src/Drivers/NucMeshDriver.C

      include/Geometry/ShapeData.H
      include/Geometry/GeoManager.H
      src/Geometry/GeoManager.C

      include/NucMesh/NucMeshGeo.H
      src/NucMesh/NucMeshGeo.C
      include/NucMesh/NucMeshShapeData.H
      src/NucMesh/NucMeshShapeData.C
      include/NucMesh/ShapeBase.H
      src/NucMesh/ShapeBase.C
      include/NucMesh/CirclesAndPolys.H
      src/NucMesh/CirclesAndPolys.C
      include/NucMesh/ShapesArray.H
      src/NucMesh/ShapesArray.C
      include/NucMesh/PolarArray.H
      src/NucMesh/PolarArray.C
      include/NucMesh/HexagonalArray.H
      src/NucMesh/HexagonalArray.C
      include/NucMesh/RectangularArray.H
      src/NucMesh/RectangularArray.C

      src/Mesh/smeshUtils.H
      include/Mesh/smeshGeoMesh.H
      src/Mesh/smeshGeoMesh.C

      include/Services/NucMeshSrv.H
      src/Services/NucMeshSrv.C
          )
endif()


if(ENABLE_OMEGAH_CUDA)
  set_source_files_properties(src/Mesh/oshGeoMesh.C src/Services/omegahRefineSrv.C
      PROPERTIES LANGUAGE CUDA)
endif()

# Building libraries ###########################################################

# Build Nemosys library
# Public headers
add_library(Nemosys ${NEMOSYS_SRCS})
target_include_directories(Nemosys
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nemosys>
)
# Private headers
target_include_directories(Nemosys PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>)

NEM_link_3p_targets(Nemosys PUBLIC ANN jsoncons omega_h eigen)
# VTK targets are irregular and can't be parsed cleanly.
target_link_libraries(Nemosys PUBLIC ${VTK_LIBRARIES})
target_compile_definitions(Nemosys PUBLIC ${VTK_DEFINITIONS})
target_link_libraries(Nemosys
    PUBLIC
        ${GMSH_LIB}
        ${EXO_LIB}
)
NEM_include_3p_directories(Nemosys PUBLIC
        ${GMSH_INCPATH}
        ${EXO_INCPATH}
)

if(WIN32)
  target_link_libraries(Nemosys PUBLIC netcdf)
endif()

if(ENABLE_GMSH)
  NEM_link_3p_targets(Nemosys PUBLIC Gmsh::Gmsh MAdLib)
endif()


if(HAVE_CGNS)
  target_link_libraries(Nemosys PUBLIC ${CGNS_LIB})
  NEM_include_3p_directories(Nemosys PUBLIC ${CGNS_INCPATH})
  target_compile_definitions(Nemosys PRIVATE HAVE_CGNS)
endif()

if(HAVE_HDF5)
  target_compile_definitions(Nemosys PRIVATE HAVE_HDF5)
  if(HAVE_HDF5_CONFIG)
    if(BUILD_SHARED_LIBS)
      NEM_link_3p_targets(Nemosys
          PUBLIC
              hdf5::hdf5_cpp-shared
              hdf5::hdf5-shared
      )
    else()
      NEM_link_3p_targets(Nemosys
          PUBLIC
              hdf5::hdf5_cpp-static
              hdf5::hdf5-static
      )
    endif()
  else()
    target_link_libraries(Nemosys PUBLIC ${HDF5_LIBRARIES})
    NEM_include_3p_directories(Nemosys PUBLIC ${HDF5_INCLUDE_DIRS})
  endif()
endif()

if(ENABLE_OPENCASCADE)
  target_link_libraries(Nemosys PUBLIC TKG3d TKBRep)
  target_link_libraries(Nemosys PRIVATE TKTopAlgo TKGeomBase TKMath TKBO)
  # OpenCASCADE doesn't set INTERFACE_INCLUDE_DIRS
  NEM_include_3p_directories(Nemosys PUBLIC ${OpenCASCADE_INCLUDE_DIR})
  NEM_link_3p_targets(Nemosys PUBLIC SMESH StdMeshers SMDS)
  NEM_link_3p_targets(Nemosys PRIVATE Controls)
  if(ENABLE_NETGEN)
    NEM_link_3p_targets(Nemosys PUBLIC NETGENPlugin)
  endif()
endif()


# TODO: Appears unused.
#if(NOT ${BUILD_SHARED_LIBS})
#  add_definitions(-DSTATIC_LINK)
#endif()

if(ENABLE_MPI)
  target_link_libraries(Nemosys PUBLIC ${MPI_LIBRARIES})
endif()

if(ENABLE_OPENMP)
  target_compile_definitions(Nemosys PUBLIC HAVE_OPENMP)
  target_link_libraries(Nemosys PUBLIC OpenMP::OpenMP_CXX)
endif()

if(ENABLE_EPIC)
  target_compile_definitions(Nemosys PUBLIC HAVE_EPIC)
  NEM_link_3p_targets(Nemosys PRIVATE kmeans_lib)

  NEM_link_3p_targets(Nemosys PUBLIC quickhull)

  NEM_link_3p_targets(Nemosys PRIVATE fast-cpp-csv-parser)
endif()

if(HAVE_METIS)
  target_link_libraries(Nemosys PUBLIC ${METIS_LIB})
  # Brought in by Gmsh already! Should avoid the clash somehow.
  target_compile_definitions(Nemosys PRIVATE HAVE_METIS)
  NEM_include_3p_directories(Nemosys PUBLIC ${METIS_INCPATH})
endif()

if(ENABLE_NETGEN)
  if(HAVE_NETGEN_OLD)
    target_link_libraries(Nemosys PUBLIC ${NETGEN_LIB})
  else()
    target_link_libraries(Nemosys PUBLIC nglib)
  endif()
  target_compile_definitions(Nemosys PUBLIC HAVE_NGEN)
  NEM_include_3p_directories(Nemosys PUBLIC ${NETGEN_INCLUDE_DIRS})
endif()

if(ENABLE_TEMPLATE_MESH)
  target_compile_definitions(Nemosys PUBLIC HAVE_TEMPLATE_MESH)
endif()

if(ENABLE_CFMSH)
  target_link_libraries(Nemosys
      PUBLIC foamutils # Link foamutils lib
      PUBLIC ${OPNF_LIBRARIES})
  NEM_include_3p_directories(Nemosys PUBLIC ${OPNF_INCLUDE_DIRS})
  NEM_include_3p_directories(Nemosys PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/contrib/cfmesh_headers>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cfmesh>
  )

  target_compile_definitions(Nemosys PUBLIC
      HAVE_CFMSH
      ${OPNF_COMPILE_DEFINITIONS}
  )

  string(FIND $ENV{WM_PROJECT_VERSION} "v2006" OF2006)
  string(FIND $ENV{WM_PROJECT_VERSION} "v2012" OF2012)
  if(OF2006 GREATER -1)
    set(OF_VERSION v2006)
    target_compile_definitions(Nemosys PUBLIC HAVE_OF2006)
  endif()
  if(OF2012 GREATER -1)
    set(OF_VERSION v2012)
    target_compile_definitions(Nemosys PUBLIC HAVE_OF2012)
  endif()

  target_link_libraries(Nemosys PRIVATE Boost::filesystem)
endif()

if(ENABLE_MLAMR)
  #set(CMAKE_CXX_STANDARD 14)
  target_compile_definitions(Nemosys PUBLIC
      MLAMR)
  find_package(frugally-deep CONFIG REQUIRED)
  message(STATUS "frugally-deep include ${PACKAGE_PREFIX_DIR}/include")
  target_include_directories(Nemosys SYSTEM
      PUBLIC ${PACKAGE_PREFIX_DIR}/include)
endif()

if(ENABLE_CONSRV_VOLUME_TRANSFER)
  target_link_libraries(Nemosys PUBLIC supermesh)
  target_compile_definitions(Nemosys PUBLIC HAVE_SUPERMESH)
endif()

if(ENABLE_CONSRV_SURFACE_TRANSFER)
    target_link_libraries(Nemosys
        PUBLIC
        IMPACT::SITCOM
        IMPACT::SurfX
        IMPACT::SimOUT)
  target_compile_definitions(Nemosys PUBLIC HAVE_IMPACT)
endif()

if(ENABLE_OPENMP)
  target_link_libraries(Nemosys
      PUBLIC
      OpenMP::OpenMP_CXX)
endif()

if(ENABLE_GMSH)
  target_compile_definitions(Nemosys PRIVATE HAVE_GMSH)
endif()

if(ENABLE_OPENCASCADE)
  target_compile_definitions(Nemosys PRIVATE HAVE_OCC)
endif()


# Set warning level ############################################################

get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if(languages MATCHES CUDA)
  set(conditional_cuda_flag "$<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=>")
else()
  set(conditional_cuda_flag "")
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  target_compile_options(Nemosys PRIVATE "${conditional_cuda_flag}-Wpedantic")
  target_compile_options(Nemosys PRIVATE "${conditional_cuda_flag}-Wall"
          "${conditional_cuda_flag}-Wno-sign-compare")
  target_compile_options(Nemosys PRIVATE "${conditional_cuda_flag}-Wextra"
          "${conditional_cuda_flag}-Wno-unused-parameter")
  target_compile_options(Nemosys PRIVATE "${conditional_cuda_flag}-Woverloaded-virtual"
          "${conditional_cuda_flag}-Wnon-virtual-dtor")
  target_compile_options(Nemosys PRIVATE "${conditional_cuda_flag}-Wstrict-null-sentinel"
          "${conditional_cuda_flag}-Wredundant-decls")

  target_compile_options(Nemosys PUBLIC "${conditional_cuda_flag}-Wno-deprecated-declarations")
#  target_compile_definitions(Nemosys PUBLIC
#      NEMOSYS_DEPRECATED_EXPORT=NEMOSYS_EXPORT
#      NEMOSYS_DEPRECATED_NO_EXPORT=NEMOSYS_NO_EXPORT)

  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)
    target_compile_options(Nemosys PRIVATE "${conditional_cuda_flag}-Wfloat-conversion")
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)

  # Stretch goals
  #target_compile_options(Nemosys PUBLIC -Wshadow)
  #target_compile_options(Nemosys PUBLIC -Wfloat-equal)
  #target_compile_options(Nemosys PUBLIC -Wold-style-cast)
  #target_compile_options(Nemosys PUBLIC -Wsign-compare -Wunused-parameter)

  # Maybe
  #target_compile_options(Nemosys PUBLIC -Wmissing-declarations)
  #target_compile_options(Nemosys PUBLIC -Wsign-promo -Wsuggest-override
  #    -Wsuggest-final-types -Wsuggest-final-methods
  #    -Wuseless-cast -Wzero-as-null-pointer-constant -Wpadded)

  # Requires high gcc
  #target_compile_options(Nemosys PUBLIC -Wredundant-tags -Wmismatched-tags
  #    -Wextra-semi -Warith-conversion
  #    -Wmultiple-inheritance -Wvirtual-inheritance)
elseif(MSVC)
  #target_compile_options(Nemosys PUBLIC /permissive-)
  if(POLICY CMP0092)
    # CMP0092 introduced in CMake 3.15. Removes /W3 as default warning flag for
    # MSVC-based compilers.
    cmake_policy(SET CMP0092 NEW)
    target_compile_options(Nemosys PUBLIC /W1)
  else()
    string(REGEX REPLACE "/W3" "/W1" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  endif()
endif()

# Generate export header #######################################################
include(GenerateExportHeader)
generate_export_header(Nemosys)
target_include_directories(Nemosys PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

# Set soname version ###########################################################

set_target_properties(Nemosys PROPERTIES VERSION ${NEMOSYS_VERSION}
                                         SOVERSION ${NEMOSYS_MAJOR_VERSION})

# Install libraries ############################################################

# Headers
# Only install public headers, not those in src
foreach(source_file ${NEMOSYS_SRCS})
  if(${source_file} MATCHES "^include/.*\\.H$")
    install(DIRECTORY include/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nemosys
            FILES_MATCHING PATTERN ${source_file})
  endif()
endforeach()

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/nemosys_export.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nemosys)

# cfmesh
if(ENABLE_CFMSH)
  install(FILES ${CFMSH_INC_FILES}
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cfmesh)
endif()

# Libraries
install(TARGETS Nemosys
    EXPORT Nemosys
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Export Nemosys as a CMake package
# Make sure to add "find_dependency(<package>)" for packages that contain targets that appear in
# INTERFACE_LINK_LIBRARIES for the Nemosys target (from either target_link_libraries(Nemosys PUBLIC|INTERFACE ...)
# or NEM_link_3p_targets(Nemosys PUBLIC|INTERFACE ...))
include(CMakePackageConfigHelpers)
set(NEM_CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/Nemosys)
configure_package_config_file(NemosysConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/NemosysConfig.cmake
    INSTALL_DESTINATION ${NEM_CONFIG_INSTALL_DIR})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NemosysConfig.cmake
    DESTINATION ${NEM_CONFIG_INSTALL_DIR})
install(EXPORT Nemosys
    NAMESPACE Nemosys::
    DESTINATION ${NEM_CONFIG_INSTALL_DIR})

# Setting utility files ########################################################

# Setting nemosysRun source
set(NEMOSYSRUN_SRCS utils/Nemosys.C)

# TODO: clean up old utilities
if(ENABLE_BUILD_UTILS)
  if(HAVE_CGNS AND ENABLE_GMSH)
    # Setting stitch cgns utility source
    set(ROCSTITCHMESH_SRCS utils/rocStitchMesh.C)
  endif()
endif()

# Building utilities ###########################################################

# Build nemosysRun executable
add_executable(nemosysRun ${NEMOSYSRUN_SRCS})
target_link_libraries(nemosysRun Nemosys)

if(ENABLE_BUILD_UTILS)
  if(HAVE_CGNS AND ENABLE_GMSH)
    # Building stitch cgns utility
    add_executable(rocStitchMesh ${ROCSTITCHMESH_SRCS})
    target_link_libraries(rocStitchMesh Nemosys)
    if(HAVE_GLOB_H)
      target_compile_definitions(rocStitchMesh PRIVATE HAVE_GLOB_H)
    else()
      target_link_libraries(rocStitchMesh Boost::filesystem)
    endif()
  endif()
endif()

# Install utilities ############################################################

install(TARGETS nemosysRun
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

if(ENABLE_BUILD_UTILS)
  if(HAVE_CGNS AND ENABLE_GMSH)
    install(TARGETS rocStitchMesh
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()
endif()

# Building Python wrapper ######################################################

if(ENABLE_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

# Building documentation #######################################################
# TODO: revamp the documentation

option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" OFF)
if(BUILD_DOCUMENTATION)
  find_package(Doxygen)
  if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
  endif()
  set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/Nemosys.cfg.in)
  set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Nemosys.cfg)
  set(documentation ${CMAKE_CURRENT_BINARY_DIR}/docs)
  configure_file(${doxyfile_in} ${doxyfile} @ONLY)
  add_custom_command(OUTPUT ${documentation}
      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
      MAIN_DEPENDENCY ${doxyfile} ${doxyfile_in}
      DEPENDS Nemosys ${doxy_extra_files}
      COMMENT "Generating documentation")
  add_custom_target(docs DEPENDS ${documentation})
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION ${CMAKE_INSTALL_DATAROOTDIR})
endif()

# Perform testing ##############################################################

if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(testing)
endif(ENABLE_TESTING)

# Packaging ####################################################################

if(WIN32)
  install(CODE "include(BundleUtilities)
  fixup_bundle(\"\${CMAKE_INSTALL_PREFIX}/bin/nemosysRun.exe\" \"\" \"${LIBRARY_DIRS}\")")
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION "bin")
  set(CMAKE_INSTALL_UCRT_LIBRARIES "TRUE")
  include(InstallRequiredSystemLibraries)
endif()

set(CPACK_PACKAGE_NAME "NEMoSys")
set(CPACK_PACKAGE_VENDOR "${NEMOSYS_AUTHOR}")
set(CPACK_PACKAGE_VERSION "${NEMOSYS_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${NEMOSYS_MAJOR_VERSION}")
set(CPACK_PACKAGE_VERSION_MINOR "${NEMOSYS_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH "${NEMOSYS_PATCH_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION "${NEMOSYS_DESCRIPTION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${NEMOSYS_DESCRIPTION}")
set(CPACK_PACKAGE_HOMEPAGE_URL "${NEMOSYS_HOMEPAGE_URL}")
set(CPACK_PACKAGE_CONTACT "${NEMOSYS_AUTHOR_EMAIL}")
#set(CPACK_PACKAGE_INSTALL_DIRECTORY "") # Use default location.
#set(CPACK_PACKAGE_ICON "") # image used by GUI installers
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")

if(WIN32)
  set(CPACK_GENERATOR "WIX")
  set(CPACK_WIX_UPGRADE_GUID "EE40DAF2-CA3C-4013-AABE-D9485EF7E692")
  set(CPACK_WIX_PRODUCT_ICON "${CMAKE_SOURCE_DIR}/nemosys_logo_color.png")
  set(CPACK_WIX_CMAKE_PACKAGE_REGISTRY "Nemosys")
else()
  if(CMAKE_VERSION VERSION_LESS 3.10)
    message(WARNING "CMake 3.10 or above is required to build a DEB package.")
  else()
    set(CPACK_GENERATOR DEB)
    #set(CPACK_DEBIAN_PACKAGE_NAME "")
    #set(CPACK_DEBIAN_PACKAGE_VERSION "")
    #set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "")
    if(ENABLE_CFMSH)
      set(CPACK_DEBIAN_PACKAGE_DEPENDS "netgen, openfoam${OF_VERSION}")
    else()
      set(CPACK_DEBIAN_PACKAGE_DEPENDS "netgen")
    endif()
    #set(CPACK_DEBIAN_PACKAGE_MAINTAINER "")
    #set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "")
    set(CPACK_DEBIAN_PACKAGE_SECTION "science")
    #set(CPACK_DEBIAN_COMPRESSION_TYPE "gzip")
    #set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
    set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "${CPACK_PACKAGE_HOMEPAGE_URL}")
    set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
    set(CPACK_DEBIAN_PACKAGE_DEBUG ON)  # For testing only
    #set(CPACK_DEBIAN_PACKAGE_PREDEPENDS "")
    #set(CPACK_DEBIAN_PACKAGE_ENHANCES "")
    #set(CPACK_DEBIAN_PACKAGE_BREAKS "")
    #set(CPACK_DEBIAN_PACKAGE_CONFLICTS "")
    #set(CPACK_DEBIAN_PACKAGE_PROVIDES "")
    #set(CPACK_DEBIAN_PACKAGE_REPLACES "")
    #set(CPACK_DEBIAN_PACKAGE_RECOMMENDS "")
    #set(CPACK_DEBIAN_PACKAGE_SUGGESTS "")
    #set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "")

    set(CPACK_DEBIAN_PACKAGE_CONTROL_STRICT_PERMISSION TRUE)  # CMake 3.4
    #set(CPACK_DEBIAN_PACKAGE_SOURCE "")  # CMake 3.5
    #set(CPACK_DEB_COMPONENT_INSTALL OFF)  # CMake 3.6
    set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")  # CMake 3.6
    set(CPACK_DEBIAN_PACKAGE_RELEASE "1")  # CMake 3.6
    set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)  # CMake 3.6
    #set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS_POLICY "=")  # CMake 3.6
    # In case of long file names and CMake < 3.14, pax archive type can cause problems
    set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar")  # CMake 3.7
    #set(CPACK_DEBIAN_PACKAGE_EPOCH "")  # Cmake 3.10

    set(CPACK_PACKAGE_RELOCATABLE FALSE)
    set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST
      /etc/ld.so.conf.d
      /etc
      /usr/local/lib64
      /usr/local
      /usr/lib/cmake
      /usr
      /usr/lib
      /usr/bin
      /usr/include
      /usr/lib64
      /usr/lib64/cmake
      /usr/lib64/pkgconfig)
  endif()
endif()
# Lists exe and text label used for Start Menu shortcuts
#set(CPACK_PACKAGE_EXECUTABLES "")
# Must have matching Start Menu shortcut to create a Desktop link
#set(CPACK_CREATE_DESKTOP_LINKS "")

# WIN32 NOTE: - Run CPack from a MSVC environment as dumpbin is required
#               for the fixup_bundle search of DLL dependencies.
#             - The system libraries pulled in are release DLLs.
#               Do NOT package in Debug! THAT APPLIES TO DEPENDENCIES TOO!!
include(CPack)
