
# TODO: Make this hierarchical
cmake_minimum_required(VERSION 3.1)

# Set project metadata #########################################################

set(NEMOSYS_MAJOR_VERSION 0)
set(NEMOSYS_MINOR_VERSION 48)
set(NEMOSYS_PATCH_VERSION 0)
set(NEMOSYS_EXTRA_VERSION "" CACHE STRING "NEMoSys extra version string")
set(NEMOSYS_VERSION "${NEMOSYS_MAJOR_VERSION}.${NEMOSYS_MINOR_VERSION}")
set(NEMOSYS_VERSION
    "${NEMOSYS_VERSION}.${NEMOSYS_PATCH_VERSION}${NEMOSYS_EXTRA_VERSION}")
set(NEMOSYS_SHORT_LICENSE "GNU General Public License v2 or later (GPLv2+)")
set(NEMOSYS_HOMEPAGE_URL "http://git.illinois.rocstar/Nemosys/Nemosys")
set(NEMOSYS_AUTHOR
    "Masoud Safdari (lead)"
    "Alessandro Gondolo"
    "Sachin Natesh"
    "Akash Patel"
    "Seth Pemberton"
    "Pavle Simonovic"
    "Taiyo Wilson")
set(NEMOSYS_AUTHOR_EMAIL "msafdari@illinoisrocstar.com")
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

# TODO: move out of main
# macro definitions

# This variable controls the default value of the options, which are normally
# set to ON (useful if you want to configure a minimal version of NEMoSys: e.g.,
#  "cmake -DDEFAULT=0 -DENABLE_DTK=1")
set(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")
# creating code compilation database
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI
      "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

list(APPEND CMAKE_BUILD_TYPE_OPTIONS
    "None" "Debug" "Release" "MinSizeRel" "RelWithDebInfo")

macro(NEM_get_imported_location OUT_VAR LIB_NAME OUT_VAR_PATH)
  foreach(buildtype IN LISTS CMAKE_BUILD_TYPE_OPTIONS)
    string(TOUPPER ${buildtype} BUILD_TYPE)
    if(WIN32)
      get_target_property(OUT_LIB ${LIB_NAME} IMPORTED_IMPLIB_${BUILD_TYPE})
    else()
      get_target_property(OUT_LIB ${LIB_NAME} IMPORTED_LOCATION_${BUILD_TYPE})
    endif()
    if(OUT_LIB)
      list(APPEND ${OUT_VAR} ${OUT_LIB})
      get_filename_component(LIBRARY_PATH ${OUT_LIB} PATH)
      list(APPEND ${OUT_VAR_PATH} ${LIBRARY_PATH})
    endif()
  endforeach()
endmacro()

# set to ${DEFAULT} to use by DEFAULT
opt(MPI "Enable MPI support" OFF)
opt(TESTING "Enable testing" ${DEFAULT})
opt(BUILD_UTILS "Build utilities" OFF)
opt(PYTHON_BINDINGS "Enable Python bindings" OFF)
opt(CFMSH "Enable cfMesh Meshing engine" OFF)
opt(CGNS "Enable CGNS extensions (Used for Rocstar utilities)" OFF)
opt(DTK "Enable DTK extensions" OFF)
opt(EXODUS "Enable EXODUS II extensions" OFF)
opt(EPIC "Enable EPIC preprocessor" OFF)
opt(HDF5 "ENABLE HDF5 extensions (Used by some CGNS and Netgen builds)" OFF)
opt(METIS "Enable METIS partitioner" ${DEFAULT})
opt(NETGEN "Enable Netgen meshing engine" ${DEFAULT})
opt(OPENCASCADE "Enable OpenCASCADE kernel support (Used by some GMSH builds)" ${DEFAULT})
opt(SIMMETRIX "Enable Simmetrix Meshing engine (Used for Rocstar remeshing)" OFF)
opt(TEMPLATE_MESH "Enable meshing templates" ON)
opt(CONSRV_SURFACE_TRANSFER "Enable Conservative Surface Transfer" OFF)
opt(CONSRV_VOLUME_TRANSFER "Enable Conservative Volume Transfer" OFF)

# Check for option interdependencies ###########################################

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
endif()

# Find VTK
# v7.1.0 required for node set and side set name support in vtkModelMetadata
find_package(VTK 7.1.0 REQUIRED)
# UseVTK.cmake file does not correctly set VTK headers as external.
#include(${VTK_USE_FILE})
set(VTK_LIB_DIR ${vtksys_RUNTIME_LIBRARY_DIRS})
list(APPEND LIBRARY_DIRS ${VTK_LIB_DIR})
list(APPEND WRAPPER_INCLUDE_DIRS ${VTK_INCLUDE_DIRS})
message(STATUS "vtk include path: ${VTK_INCLUDE_DIRS}")
message(STATUS "vtk libs ${VTK_LIB_DIR}")

# Finding vtkCommonCore library to pass to Python for linking
get_target_property(VTK_COMMON_CORE_ILL vtkCommonCore INTERFACE_LINK_LIBRARIES)
if(VTK_COMMON_CORE_ILL)
  # For VTK 7
  NEM_get_imported_location(VTK_COMMON_CORE_LIB ${VTK_COMMON_CORE_ILL} VTK_COMMON_CORE_LIB_PATH)
else()
  # For VTK 6 and 8
  NEM_get_imported_location(VTK_COMMON_CORE_LIB vtkCommonCore VTK_COMMON_CORE_LIB_PATH)
endif()
get_filename_component(VTK_COMMON_CORE_LIB ${VTK_COMMON_CORE_LIB} NAME)
if(WIN32) # Remove the .lib suffix
  string(REGEX REPLACE "\\.lib$" "" VTK_COMMON_CORE_LIB ${VTK_COMMON_CORE_LIB})
endif()
if(CMAKE_COMPILER_IS_GNUCXX) # Add the : in front to force exact match by ld
  set(VTK_COMMON_CORE_LIB ":${VTK_COMMON_CORE_LIB}")
endif()
message(STATUS "vtkCommonCore library ${VTK_COMMON_CORE_LIB}")

# Find Gmsh
find_library(GMSH_LIB gmsh HINTS lib)
find_file(GMSH_HDR gmsh.h HINTS include)
get_filename_component(GMSH_INCPATH ${GMSH_HDR} PATH)
get_filename_component(GMSH_LIB_PATH ${GMSH_LIB} PATH)
list(APPEND LIBRARY_DIRS ${GMSH_LIB_PATH})
list(APPEND WRAPPER_INCLUDE_DIRS ${GMSH_INCPATH})
message(STATUS "gmsh include path ${GMSH_INCPATH}")
message(STATUS "gmsh libs ${GMSH_LIB}")

# Find OpenCascade (used by Gmsh if build without ENABLE_OCC_STATIC)
if(ENABLE_OPENCASCADE)
  find_package(OpenCASCADE REQUIRED)
  list(APPEND LIBRARY_DIRS ${OpenCASCADE_LIBRARY_DIR})
  list(APPEND WRAPPER_INCLUDE_DIRS ${OpenCASCADE_INCLUDE_DIR})
  message(STATUS "OpenCASCADE include path ${OpenCASCADE_INCLUDE_DIR}")
  message(STATUS "OpenCASCADE libs ${OpenCASCADE_LIBRARY_DIR}")

  # CMake 3.13 required for target_link_directories
  link_directories(${OpenCASCADE_LIBRARY_DIR})
  # /CMake 3.13

  set(HAVE_OPENCASCADE TRUE)
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
  list(APPEND WRAPPER_INCLUDE_DIRS ${CGNS_INCPATH})
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
    list(APPEND WRAPPER_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
    message(STATUS "hdf5 include path ${HDF5_INCLUDE_DIRS}")
    message(STATUS "hdf5 libs ${HDF5_LIBRARIES}")
    set(HAVE_HDF5 TRUE)
    set(HAVE_HDF5_CONFIG TRUE)
  else()
    find_package(HDF5 REQUIRED COMPONENTS CXX)
    list(GET ${HDF5_LIBRARIES} 0 HDF5_LIBRARY)
    get_filename_component(HDF5_LIB_PATH ${HDF5_LIBRARY} PATH)
    list(APPEND LIBRARY_DIRS ${HDF5_LIB_PATH})
    list(APPEND WRAPPER_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
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
  list(APPEND WRAPPER_INCLUDE_DIRS ${METIS_INCPATH})
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
  list(APPEND WRAPPER_INCLUDE_DIRS ${NETGEN_INCLUDE_DIRS})
  message(STATUS "netgen include path ${NETGEN_INCLUDE_DIRS}")
  message(STATUS "netgen libs ${NETGEN_LIBRARY_DIR}")

  # Find zlib
  if(WIN32)
    add_definitions(-D_AFXDLL) # Needed by Netgen headers.
    find_package(ZLIB REQUIRED) # Find zlib headers since Netgen doesn't
    include_directories(${ZLIB_INCLUDE_DIRS})
    get_filename_component(ZLIB_LIB_PATH ${ZLIB_LIBRARIES} PATH)
    list(APPEND LIBRARY_DIRS ${ZLIB_LIB_PATH} ${ZLIB_LIB_PATH}/../bin)
    list(APPEND WRAPPER_INCLUDE_DIRS ${ZLIB_INCLUDE_DIRS})
    message(STATUS "zlib include path ${ZLIB_INCLUDE_DIRS}")
    message(STATUS "zlib libs ${ZLIB_LIBRARIES}")
  endif()
endif()

# Find ExodusII
if(ENABLE_EXODUS)
  find_library(EXO_LIB NAMES exoIIv2 exodus exoIIv2c)
  find_file(EXO_HDR exodusII.h)
  get_filename_component(EXO_LIB_PATH ${EXO_LIB} PATH)
  get_filename_component(EXO_INCPATH ${EXO_HDR} PATH)
  list(APPEND LIBRARY_DIRS ${EXO_LIB_PATH})
  list(APPEND WRAPPER_INCLUDE_DIRS ${EXO_INCPATH})
  message(STATUS "exodus include path ${EXO_INCPATH}")
  message(STATUS "exodus libs ${EXO_LIB}")

  # Find NetCDF
  if(WIN32)
    find_package(netCDF REQUIRED)
    list(APPEND LIBRARY_DIRS ${netCDF_LIB_DIR} ${netCDF_LIB_DIR}/../bin)
    list(APPEND WRAPPER_INCLUDE_DIRS ${netCDF_INCLUDE_DIRS})
    message(STATUS "netcdf include path ${netCDF_INCLUDE_DIR}")
    message(STATUS "netcdf libs ${netCDF_LIB_DIR}")
  endif()
endif()

# Find Trilinos
if(ENABLE_DTK)
  find_package(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos
                              ${TRILINOS_PATH})
  if(NOT Trilinos_FOUND)
    message(FATAL_ERROR "Could not find Trilinos!")
  endif()
  message(Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
  link_directories(${Trilinos_TPL_LIBRARY_DIRS})

  # Import dtk related variables
  include(dtk_stuff.txt)
endif()

# Find Simmetrix
if(ENABLE_SIMMETRIX)
  set(SIMMETRIX_LIB_DIR CACHE PATH "Simmetrix library directory path.")
  find_library(SIMADM_LIB SimAdvMeshing HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMEXP_LIB SimExport HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMDSC_LIB SimDiscrete HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMPMH_LIB SimPartitionedMesh HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMPWR_LIB SimPartitionWrapper HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMMSH_LIB SimMeshing HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMMTL_LIB SimMeshTools HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMMDL_LIB SimModel HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMLCN_LIB SimLicense HINTS ${SIMMETRIX_LIB_DIR})
  set(SIMMETRIX_LIBS ${SIMADM_LIB} ${SIMEXP_LIB} ${SIMDSC_LIB} ${SIMPMH_LIB}
      ${SIMPWR_LIB} ${SIMMSH_LIB} ${SIMMTL_LIB} ${SIMMDL_LIB} ${SIMLCN_LIB})
  find_file(SIMMETRIX_HDR MeshSim.h HINTS ${SIMMETRIX_LIB_DIR}/../../include)
  get_filename_component(SIMMETRIX_INCPATH ${SIMMETRIX_HDR} PATH)
  get_filename_component(SIMMETRIX_LIB_PATH ${SIMMDL_LIB} PATH)
  list(APPEND WRAPPER_INCLUDE_DIRS ${SIMMETRIX_INCPATH})
  message(STATUS "Simmetrix include path ${SIMMETRIX_INCPATH}")
  message(STATUS "Simmetrix libs ${SIMMETRIX_LIBS}")
endif()

# Find boost
# Needed if there is no glob.h on the system
if(ENABLE_SIMMETRIX OR ENABLE_CFMSH OR NOT HAVE_GLOB_H)
  if(WIN32)
    set(Boost_USE_STATIC_LIBS ON) # needed ON on WIN32.
  endif()
  set(Boost_USE_MULTITHREADED ON)
  find_package(Boost 1.45.0 REQUIRED COMPONENTS filesystem regex)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    list(APPEND LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
    list(APPEND WRAPPER_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
  endif()
endif()

# Find cfMesh and OpenFOAM
if(ENABLE_CFMSH)
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

  message(STATUS "OpenFOAM include directory ${OPNF_INCLUDE_DIRS}")
  message(STATUS "OpenFOAM link directory ${OPNF_LIBRARY_DIRS}")

  find_package(OpenMP REQUIRED)

  # CMake 3.13 required for target_link_options
  set(CMAKE_EXE_LINKER_FLAGS ${OpenMP_CXX_FLAGS})
  set(CMAKE_MODULE_LINKER_FLAGS ${OpenMP_CXX_FLAGS})
  set(CMAKE_SHARED_LINKER_FLAGS ${OpenMP_CXX_FLAGS})
  set(CMAKE_STATIC_LINKER_FLAGS ${OpenMP_CXX_FLAGS})
  # /CMake 3.13
endif()

if(ENABLE_CONSRV_SURFACE_TRANSFER)
  find_package(IMPACT REQUIRED)
endif()

# For some reason, older versions of cmake fail to find the local gmsh include
# files when compiling MAdLib header references.
# TODO: Is this still true in Gmsh 4.2+ and CMake 2.8+?
#if(EXISTS ${GMSH_INCPATH})
#  include_directories(AFTER ${GMSH_INCPATH}/../)
#endif()

# Setting compile and link options #############################################

# Adding flags to Debug
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# Compiled objects output paths
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Find contrib libraries #######################################################

# Find MAdLib and ANN
add_subdirectory(contrib/madlib-1.3.0) # includes ANN
list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0)
list(APPEND LIBRARY_DIRS ${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0/Contrib/ANN)
get_target_property(MADLIB_INCPATH MAdLib INTERFACE_INCLUDE_DIRECTORIES)
list(APPEND WRAPPER_INCLUDE_DIRS ${MADLIB_INCPATH})
get_target_property(ANN_INCPATH ANN INTERFACE_INCLUDE_DIRECTORIES)
list(APPEND WRAPPER_INCLUDE_DIRS ${ANN_INCPATH})

# Find jsoncons
add_subdirectory(contrib/jsoncons-0.127.0)
get_target_property(JSONCONSRV_INCPATH jsoncons INTERFACE_INCLUDE_DIRECTORIES)
list(APPEND WRAPPER_INCLUDE_DIRS ${JSONCONSRV_INCPATH})

if(ENABLE_CFMSH)
  # Find cfMesh
  add_subdirectory(contrib/cfmesh-1.1.2ir1)
  # TODO: Add LIBRARY_DIRS and WRAPPER_INCLUDE_DIRS
endif()

if(ENABLE_CONSRV_VOLUME_TRANSFER)
  add_subdirectory(contrib/libsupermesh-1.0.1)
  # TODO: Add LIBRARY_DIRS and WRAPPER_INCLUDE_DIRS
endif()

# Setting library files ########################################################
# TODO: clean up the unnecessary modules such as R8

# Setting include files
file(GLOB H_FILES include/*.H)
file(GLOB h_FILES include/*.h)
file(GLOB hpp_FILES include/*.hpp)
file(GLOB cpp_FILES include/*.cpp)
set(Nemosys_INC_FILES ${H_FILES} ${h_FILES} ${hpp_FILES} ${cpp_FILES})

# Setting r8 source
set(R8_SRCS
    src/Math/r8lib.cpp
)

# Setting interp source
set(INTERP_SRCS
    src/Geometry/spheres.C

    src/Math/baseInterp.C
    src/Math/rbf_interp_nd.cpp
    src/Math/rbfInterp.C
)

# Setting Nemosys source
# TODO: Commented files are orphaned! Check if to delete or reinclude.
set(NEMOSYS_SRCS
    src/Drivers/NemDriver.C
    src/Drivers/ConversionDriver.C
    src/Drivers/InputGenDriver.C
    src/Drivers/MeshGenDriver.C
    src/Drivers/MeshQualityDriver.C
    src/Drivers/NucMeshDriver.C
    src/Drivers/RefineDriver.C
    src/Drivers/TransferDriver.C
    src/Drivers/PackMeshDriver.C

    src/Geometry/convexContainer.C
    src/Geometry/point.C
    src/Geometry/qhQuickHull.C
    src/Geometry/rocPack.C
    src/Geometry/rocPackShape.C
    src/Geometry/hmxShape.C
    src/Geometry/petnShape.C
    src/Geometry/icosidodecahedronShape.C

    src/InputGeneration/inputGen.C

    src/Integration/Cubature.C

    src/Math/kmeans.C

    src/Mesh/meshBase.C
    src/Mesh/cobalt.C
    src/Mesh/gmshMesh.C
    src/Mesh/patran.C
    src/Mesh/pntMesh.C
    src/Mesh/vtkMesh.C

    src/MeshGeneration/meshGen.C
    src/MeshGeneration/netgenGen.C
    src/MeshGeneration/netgenParams.C

    src/MeshPartitioning/meshPartitioner.C

    src/MeshQuality/MeshQuality.C

    src/NucMesh/circles.C
    src/NucMesh/polygon.C
    src/NucMesh/shape.C

    src/PatchRecovery/orthoPoly1D.C
    src/PatchRecovery/orthoPoly3D.C
    src/PatchRecovery/patchRecovery.C
    src/PatchRecovery/polyApprox.C

    src/Refinement/Refine.C

    src/SizeFieldGeneration/SizeFieldGen.C
    src/SizeFieldGeneration/GradSizeField.C
    src/SizeFieldGeneration/ValSizeField.C
    src/SizeFieldGeneration/Z2ErrorSizeField.C

    src/SolutionVerification/OrderOfAccuracy.C
#    src/SolutionVerification/RichardsonExtrapolation.C

    src/Transfer/TransferBase.C
    src/Transfer/FETransfer.C
#   src/Transfer/ConservativeVolumeTransfer.C

    src/StlToVtk.C
    src/IO/vtkAnalyzer.C
)

# Optional source files
if(HAVE_CGNS)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Drivers/RocPartCommGenDriver.C
      # src/Drivers/RocPrepDriver.C
      # src/Drivers/RocRestartDriver.C

      src/MeshPartitioning/meshStitcher.C

      src/IO/cgnsAnalyzer.C
      src/IO/cgnsWriter.C
      src/gridTransfer.C
      src/IO/rocstarCgns.C
  )
endif()

if(HAVE_HDF5)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Drivers/ProteusDriver.C

      src/IO/hdf5Reader.C
      src/IO/proteusHdf5.C
  )
endif()

if(ENABLE_EXODUS)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Mesh/exoMesh.C

      src/MeshOperation/meshSrch.C
  )
endif()

if(ENABLE_EPIC)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/InputGeneration/ep16Post.C
      src/InputGeneration/ep16Prep.C
  )
endif()

if(ENABLE_SIMMETRIX)
  if(NOT ENABLE_CGNS)
      MESSAGE(FATAL_ERROR "If enabling simmetrix, you MUST enable cgns as well.")
  endif(NOT ENABLE_CGNS)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Drivers/RemeshDriver.C

      src/MeshGeneration/simmetrixGen.C
      src/MeshGeneration/simmetrixParams.C
  )
endif()

if (ENABLE_TEMPLATE_MESH)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Drivers/TemplateMeshDriver.C
  )
endif()


if(ENABLE_CFMSH)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Mesh/foamMesh.C

      src/MeshGeneration/blockMeshGen.C
      src/MeshGeneration/blockMeshParams.C
      src/MeshGeneration/cfmeshGen.C
      src/MeshGeneration/cfmeshParams.C
      src/MeshGeneration/snappymeshGen.C
      src/MeshGeneration/snappymeshParams.C

      src/MeshManipulationFoam/MeshManipulationFoam.C
      src/MeshManipulationFoam/MeshManipulationFoamParams.C
      src/MeshManipulationFoam/mergePolyMesh.C

      src/MeshQuality/cfmshQualityParams.C
  )
endif()

if(ENABLE_CONSRV_VOLUME_TRANSFER)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
    src/Transfer/ConservativeVolumeTransfer.C
  )
endif()

if(ENABLE_CONSRV_SURFACE_TRANSFER)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Transfer/ConservativeSurfaceTransfer.C
  )
endif()

# Building libraries ###########################################################

# Build r8 and rbf library
add_library(r8 ${R8_SRCS})
target_include_directories(r8
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nemosys>
)

# Build interp library
add_library(interp ${INTERP_SRCS})
target_link_libraries(interp
    PRIVATE
        r8
        ANN
)
target_include_directories(interp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nemosys>
)

# Build Nemosys library
add_library(Nemosys ${NEMOSYS_SRCS})
# /* Work around to disable warnings coming from contrib headers.
foreach(lib jsoncons ANN MAdLib)
  get_target_property(${lib}_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
  target_include_directories(Nemosys SYSTEM PUBLIC ${${lib}_include_dirs})
endforeach()
# */
target_link_libraries(Nemosys
    PUBLIC
        interp
        MAdLib
        jsoncons
        ${VTK_LIBRARIES}
        ${GMSH_LIB}
        ${Boost_LIBRARIES}
)
target_include_directories(Nemosys
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/nemosys>
)
target_include_directories(Nemosys SYSTEM
    PUBLIC
        ${VTK_INCLUDE_DIRS}
        ${GMSH_INCPATH}
)
target_compile_definitions(Nemosys
    PUBLIC
        ${VTK_DEFINITIONS}
)

if(HAVE_CGNS)
  target_link_libraries(Nemosys PUBLIC ${CGNS_LIB})
  target_include_directories(Nemosys SYSTEM PUBLIC ${CGNS_INCPATH})
  target_compile_definitions(Nemosys PRIVATE HAVE_CGNS)
endif()

if(HAVE_HDF5)
  target_compile_definitions(Nemosys PRIVATE HAVE_HDF5)
  if(HAVE_HDF5_CONFIG)
    if(BUILD_SHARED_LIBS)
      target_link_libraries(Nemosys
          PUBLIC
              hdf5::hdf5_cpp-shared
              hdf5::hdf5-shared
      )
    else()
      target_link_libraries(Nemosys
          PUBLIC
              hdf5::hdf5_cpp-static
              hdf5::hdf5-static
      )
    endif()
  else()
    target_link_libraries(Nemosys PUBLIC ${HDF5_LIBRARIES})
    target_include_directories(Nemosys SYSTEM PUBLIC ${HDF5_INCLUDE_DIR})
  endif()
endif()

if(HAVE_OPENCASCADE)
  target_link_libraries(Nemosys PUBLIC TKernel)
endif()

if(HAVE_GLOB_H)
  target_compile_definitions(Nemosys PUBLIC HAVE_GLOB_H)
endif()

# TODO: Appears unused.
#if(NOT ${BUILD_SHARED_LIBS})
#  add_definitions(-DSTATIC_LINK)
#endif()

if(ENABLE_MPI)
  target_link_libraries(Nemosys PUBLIC ${MPI_LIBRARIES})
endif()

if(ENABLE_EXODUS)
  target_link_libraries(Nemosys PUBLIC ${EXO_LIB})
  target_compile_definitions(Nemosys PUBLIC HAVE_EXODUSII)
  target_include_directories(Nemosys SYSTEM PUBLIC ${EXO_INCPATH})

  if(WIN32)
    target_link_libraries(Nemosys PUBLIC netcdf)
  endif()
endif()

if(ENABLE_EPIC)
  target_compile_definitions(Nemosys PUBLIC HAVE_EPIC)
endif()

if(HAVE_METIS)
  target_link_libraries(Nemosys PUBLIC ${METIS_LIB})
  # Brought in by Gmsh already! Should avoid the clash somehow.
  target_compile_definitions(Nemosys PRIVATE HAVE_METIS)
  target_include_directories(Nemosys SYSTEM PUBLIC ${METIS_INCPATH})
endif()

if(ENABLE_NETGEN)
  if(HAVE_NETGEN_OLD)
    target_link_libraries(Nemosys PUBLIC ${NETGEN_LIB})
  else()
    target_link_libraries(Nemosys PUBLIC nglib)
  endif()
  target_compile_definitions(Nemosys PUBLIC HAVE_NGEN)
  target_include_directories(Nemosys SYSTEM PUBLIC ${NETGEN_INCLUDE_DIRS})
endif()

if(ENABLE_SIMMETRIX)
  target_link_libraries(Nemosys PUBLIC ${SIMMETRIX_LIBS})
  target_include_directories(Nemosys SYSTEM PUBLIC ${SIMMETRIX_INCPATH})
  target_compile_definitions(Nemosys PUBLIC HAVE_SIMMETRIX)
endif()

if(ENABLE_TEMPLATE_MESH)
  target_compile_definitions(Nemosys PUBLIC HAVE_TEMPLATE_MESH)
endif()

if(ENABLE_CFMSH)
  target_link_libraries(Nemosys
      # TODO: Add namespace to cfMesh to protect the symbol clash between
      #       libdynamicMesh.so in OpenFOAM and libmeshLibrary.so in cfMesh.
      #       Till then, link cfMesh before OpenFOAM.
      PUBLIC cfmesh
      PUBLIC ${OPNF_LIBRARIES})
  target_include_directories(Nemosys SYSTEM
      PUBLIC ${OPNF_INCLUDE_DIRS})
  target_compile_definitions(Nemosys PUBLIC
      HAVE_CFMSH
      ${OPNF_COMPILE_DEFINITIONS}
  )

  string(FIND $ENV{WM_PROJECT_VERSION} "4" OF4)
  if(OF4 GREATER -1)
    set(OF_VERSION 4)
    target_compile_definitions(Nemosys PUBLIC HAVE_OF4)
  endif()

  string(FIND $ENV{WM_PROJECT_VERSION} "5" OF5)
  if(OF5 GREATER -1)
    set(OF_VERSION 5)
    target_compile_definitions(Nemosys PUBLIC HAVE_OF5)
  endif()

  string(FIND $ENV{WM_PROJECT_VERSION} "6" OF6)
  if(OF6 GREATER -1)
    set(OF_VERSION 6)
    target_compile_definitions(Nemosys PUBLIC HAVE_OF6)
  endif()

  string(FIND $ENV{WM_PROJECT_VERSION} "7" OF7)
  if(OF7 GREATER -1)
    set(OF_VERSION 7)
    target_compile_definitions(Nemosys PUBLIC HAVE_OF7)
  endif()
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

# Set warning level ############################################################

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  target_compile_options(Nemosys PUBLIC -Wpedantic)
  target_compile_options(Nemosys PUBLIC -Wall -Wno-sign-compare)
  target_compile_options(Nemosys PUBLIC -Wextra -Wno-unused-parameter)
  target_compile_options(Nemosys PUBLIC -Woverloaded-virtual -Wnon-virtual-dtor)
  target_compile_options(Nemosys PUBLIC -Wstrict-null-sentinel -Wredundant-decls)
  if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)
      target_compile_options(Nemosys PUBLIC -Wfloat-conversion)
  endif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.9)

  # Stretch goals
  #target_compile_options(Nemosys PUBLIC -Wshadow)
  #target_compile_options(Nemosys PUBLIC -Wold-style-cast)
  #target_compile_options(Nemosys PUBLIC -Wsign-compare -Wunused-parameter)

  # Maybe
  #target_compile_options(Nemosys PUBLIC -Wmissing-declarations)
  #target_compile_options(Nemosys PUBLIC -Wsign-promo -Wsuggest-override
  #    -Wsuggest-final-types -Wsuggest-final-methods -Wfloat-equal
  #    -Wuseless-cast -Wzero-as-null-pointer-constant -Wpadded)

  # Requires high gcc
  #target_compile_options(Nemosys PUBLIC -Wredundant-tags -Wmismatched-tags
  #    -Wextra-semi -Warith-conversion
  #    -Wmultiple-inheritance -Wvirtual-inheritance)
elseif(MSVC)
#  target_compile_options(Nemosys PUBLIC /permissive-)
  string(REGEX REPLACE "/W3" "/W1" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
endif()

# Generate export header #######################################################
include(GenerateExportHeader)
generate_export_header(Nemosys)
generate_export_header(interp)
generate_export_header(r8)
list(APPEND ${Nemosys_INC_FILES} Nemosys_export.h)
list(APPEND ${interp_INC_FILES} interp_export.h)
list(APPEND ${r8_INC_FILES} r8_export.h)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/nemosys_export.h
          ${CMAKE_CURRENT_BINARY_DIR}/interp_export.h
          ${CMAKE_CURRENT_BINARY_DIR}/r8_export.h
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Set soname version ###########################################################

# TODO: Set r8 and interp versions properly. They should be in contrib!
set_target_properties(r8 PROPERTIES VERSION ${NEMOSYS_VERSION}
                                    SOVERSION ${NEMOSYS_MAJOR_VERSION})
set_target_properties(interp PROPERTIES VERSION ${NEMOSYS_VERSION}
                                        SOVERSION ${NEMOSYS_MAJOR_VERSION})
set_target_properties(Nemosys PROPERTIES VERSION ${NEMOSYS_VERSION}
                                         SOVERSION ${NEMOSYS_MAJOR_VERSION})

# Install libraries ############################################################

# Headers
install(FILES ${r8_INC_FILES} ${interp_INC_FILES} ${Nemosys_INC_FILES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nemosys)

# Libraries
install(TARGETS r8 interp Nemosys
    EXPORT Nemosys
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Export Nemosys as a CMake package
install(FILES NemosysConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Nemosys)
install(EXPORT Nemosys
    NAMESPACE Nemosys::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Nemosys)

# Setting utility files ########################################################

# Setting nemosysRun source
set(NEMOSYSRUN_SRCS utils/Nemosys.C)

# TODO: clean up old utilities
if(ENABLE_BUILD_UTILS)
  # Setting xmlDump utility source
  set(XMLDUMP_SRCS utils/xmlDump.C)

  if(HAVE_CGNS)
    # Setting cgns2gmsh utility source
    set(CGNS2GMSH_SRCS utils/cgns2msh.C)

    # Setting stitch cgns utility source
    set(ROCSTITCHMESH_SRCS utils/rocStitchMesh.C)

    # Setting Rocstar volume remesh utility source
    set(ROCREMESH_SRCS utils/rocRemesh.C)

    # Setting Rocstar surface remesh utility source
    set(ROCSURFREMESH_SRCS utils/rocSurfRemesh.C)

    # Setting Rocstar surface remesh utility source
    set(GRID2GRIDTRANSFER_SRCS utils/grid2gridTransfer.C)

    # Setting meshTransfer tutorial source
    set(MSHTRANSFER_SRCS tutorials/mshTransfer.C)
  endif()

  # TODO: Orphaned files! Check if used.
#  set(TUTR_SRCS tutorials/writeStructured.C)
#  set(UTIL_SRCS
#      utils/refineMesh.C
#      utils/vol2planeTransfer.C
#  )

  # Setting DTK tutorial source
  if(ENABLE_DTK)
    set(DTKTUTORIAL_SRC tutorials/dtkTutorial.C)
  endif()
endif()

# Building utilities ###########################################################

# Build nemosysRun executable
add_executable(nemosysRun ${NEMOSYSRUN_SRCS})
target_link_libraries(nemosysRun Nemosys)

if(ENABLE_BUILD_UTILS)
  # Building xmlDump utility
  add_executable(xmlDump MACOSX_BUNDLE ${XMLDUMP_SRCS})
  target_link_libraries(xmlDump Nemosys)
  if(NOT VTK_LIBRARIES)
    target_link_libraries(xmlDump vtkHybrid vtkWidgets)
  endif()

  if(HAVE_CGNS)
    # Building cgns2gmsh utility
    add_executable(cgns2msh ${CGNS2GMSH_SRCS})
    target_link_libraries(cgns2msh Nemosys)

    # Building stitch cgns utility
    add_executable(rocStitchMesh ${ROCSTITCHMESH_SRCS})
    target_link_libraries(rocStitchMesh Nemosys)

    # Building Rocstar volume remesh utility
    add_executable(rocRemesh ${ROCREMESH_SRCS})
    target_link_libraries(rocRemesh Nemosys)

    # Building Rocstar surface remesh utility
    add_executable(rocSurfRemesh ${ROCSURFREMESH_SRCS})
    target_link_libraries(rocSurfRemesh Nemosys)

    # Building Rocstar surface remesh utility
    add_executable(grid2gridTransfer ${GRID2GRIDTRANSFER_SRCS})
    target_link_libraries(grid2gridTransfer Nemosys)

    # Building meshTransfer tutorial
    add_executable(mshTransfer ${MSHTRANSFER_SRCS})
    target_link_libraries(mshTransfer Nemosys)
  endif()

  # Building DTK tutorial
  if(ENABLE_DTK)
    add_executable(dtkTutorial ${DTKTUTORIAL_SRC})
    set_target_properties(dtkTutorial PROPERTIES
        COMPILE_FLAGS ${dtk_headers}
        LINK_FLAGS ${dtk_libs})
    target_link_libraries(dtkTutorial Nemosys
        ${DTK_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    install(TARGETS dtkTutorial
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()
endif()

# Install utilities ############################################################

install(TARGETS nemosysRun
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

if(ENABLE_BUILD_UTILS)
  install(TARGETS xmlDump
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  if(HAVE_CGNS)
    install(TARGETS grid2gridTransfer
                    cgns2msh
                    rocStitchMesh
                    rocRemesh
                    rocSurfRemesh
                    mshTransfer
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()

  if(ENABLE_DTK)
    install(TARGETS dtkTutorial
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
set(CPACK_PACKAGE_HOMEPAGE_URL "")
set(CPACK_PACKAGE_CONTACT "${NEMOSYS_AUTHOR_EMAIL}")
#set(CPACK_PACKAGE_INSTALL_DIRECTORY "") # Use default location.
#set(CPACK_PACKAGE_ICON "") # image used by GUI installers
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")

if(WIN32)
  set(CPACK_GENERATOR NSIS)
else()
  if(CMAKE_VERSION VERSION_LESS 3.10)
    message(WARNING "CMake 3.10 or above is required to build a DEB package.")
  else()
    set(CPACK_GENERATOR DEB)
    #set(CPACK_DEBIAN_PACKAGE_NAME "")
    #set(CPACK_DEBIAN_PACKAGE_VERSION "")
    #set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "")
    if(ENABLE_CFMSH)
      set(CPACK_DEBIAN_PACKAGE_DEPENDS "libnglib-6.2.1804, openfoam${OF_VERSION}")
    else()
      set(CPACK_DEBIAN_PACKAGE_DEPENDS "libnglib-6.2.1804")
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
    #set(CPACK_DEBIAN_PACKAGE_RELEASE "1")  # CMake 3.6
    set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)  # CMake 3.6
    #set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS_POLICY "=")  # CMake 3.6
    #set(CPACK_DEBIAN_ARCHIVE_TYPE "gnutar")  # CMake 3.7
    #set(CPACK_DEBIAN_PACKAGE_EPOCH "")  # Cmake 3.10
  endif()
endif()
# Lists exe and text label used for Start Menu shortcuts
#set(CPACK_PACKAGE_EXECUTABLES "")
# Must have matching Start Menu shortcut to create a Desktop link
#set(CPACK_CREATE_DESKTOP_LINKS "")

# WIN32 NOTE: - Run CPack from a MSVC environment as dumpbin is required
#               for the fixup_bundle search of DLL dependecies.
#             - The system libraries pulled in are release DLLs.
#               Do NOT package in Debug! THAT APPLIES TO DEPENDENCIES TOO!!
include(CPack)
