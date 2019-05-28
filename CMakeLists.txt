cmake_minimum_required(VERSION 2.8)
project(Nemosys CXX C)

# This variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of NEMoSys: e.g.,
#  "cmake -DDEFAULT=0 -DENABLE_DTK=1")
set(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")
# creating code compilation database
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Setting options ##############################################################

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI
      "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

# set to ${DEFAULT} to use by DEFAULT
opt(MPI "Enable MPI support" OFF)
opt(DTK "Enable DTK extensions" OFF)
opt(EXODUS "Enable EXODUS II extensions" OFF)
opt(EPIC "Enable EPIC preprocessor" OFF)
opt(PYTHON_BINDINGS "Enable Python bindings" OFF)
opt(TESTING "Enable testing" ${DEFAULT})
opt(BUILD_UTILS "Build utilities" OFF)
opt(METIS "Enable Metis interface" ${DEFAULT})
opt(NETGEN "Enable Netgen interface" ${DEFAULT})
opt(SIMMETRIX "Enable Simmetrix Meshing engine" OFF)

# Set project metadata #########################################################

set(NEMOSYS_MAJOR_VERSION 0)
set(NEMOSYS_MINOR_VERSION 29)
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

# Set project defaults #########################################################

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
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
set(CMAKE_INSTALL_RPATH "$ORIGIN")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
if(WIN32)
  set(VTK_BUILD_SHARED_LIBS OFF)
endif()
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
get_target_property(VTK_LIB_DIR vtksys
    IMPORTED_LOCATION_${CMAKE_BUILD_TYPE_UPPER})
get_filename_component(VTK_LIB_DIR ${VTK_LIB_DIR} PATH)
list(APPEND LIBRARY_DIRS ${VTK_LIB_DIR})
message(STATUS "vtk include path: ${VTK_INCLUDE_DIRS}")
message(STATUS "vtk libs ${VTK_LIB_DIR}")

# Find MAdLib
add_subdirectory(contrib/madlib-1.3.0)
get_target_property(MADLIB_INCPATH MAdLib INTERFACE_INCLUDE_DIRECTORIES)

# Find Gmsh
find_library(GMSH_LIB gmsh HINTS lib)
find_file(GMSH_HDR gmsh.h HINTS include)
get_filename_component(GMSH_INCPATH ${GMSH_HDR} PATH)
get_filename_component(GMSH_LIB_PATH ${GMSH_LIB} PATH)
list(APPEND LIBRARY_DIRS ${GMSH_LIB_PATH})
message(STATUS "gmsh include path ${GMSH_INCPATH}")
message(STATUS "gmsh libs ${GMSH_LIB}")
add_definitions(-D_HAVE_GMSH_)

# Find OpenCascade (used by Gmsh)
find_package(OpenCASCADE REQUIRED)
message(STATUS "OpenCASCADE include path ${OpenCASCADE_INCLUDE_DIR}")
message(STATUS "OpenCASCADE libs ${OpenCASCADE_LIBRARY_DIR}")

# Find HDF5
if(WIN32)
  set(HDF5_USE_STATIC_LIBRARIES ON) # Needed on WIN32. TODO: Test on UNIX
endif()
find_package(HDF5 REQUIRED)
list(GET ${HDF5_LIBRARIES} 0 HDF5_LIBRARY)
get_filename_component(HDF5_LIB_PATH ${HDF5_LIBRARY} PATH)
list(APPEND LIBRARY_DIRS ${HDF5_LIB_PATH})
message(STATUS "hdf5 include path ${HDF5_INCLUDE_DIRS}")
message(STATUS "hdf5 libs ${HDF5_LIBRARIES}")

# Find CGNS
find_library(CGNS_LIB cgns HINTS ${CMAKE_INSTALL_PREFIX}/cgns/lib)
find_file(CGNS_HDR cgnslib.h HINTS ${CMAKE_INSTALL_PREFIX}/cgns/include)
get_filename_component(CGNS_INCPATH ${CGNS_HDR} PATH)
get_filename_component(CGNS_LIB_PATH ${CGNS_LIB} PATH)
list(APPEND LIBRARY_DIRS ${CGNS_LIB_PATH})
message(STATUS "cgns include path ${CGNS_INCPATH}")
message(STATUS "cgns libs ${CGNS_LIB}")

# Find METIS
if(ENABLE_METIS)
  find_library(METIS_LIB metis HINTS ${CMAKE_INSTALL_PREFIX}/metis/lib)
  find_file(METIS_HDR metis.h HINTS ${CMAKE_INSTALL_PREFIX}/metis/include)
  get_filename_component(METIS_INCPATH ${METIS_HDR} PATH)
  get_filename_component(METIS_LIB_PATH ${METIS_LIB} PATH)
  list(APPEND LIBRARY_DIRS ${METIS_LIB_PATH})
  message(STATUS "metis include path ${METIS_INCPATH}")
  message(STATUS "metis libs ${METIS_LIB}")
endif()

# Find Netgen
if(ENABLE_NETGEN)
  find_package(Netgen REQUIRED)
  list(APPEND LIBRARY_DIRS ${NETGEN_LIBRARY_DIR})
  message(STATUS "netgen include path ${NETGEN_INCLUDE_DIRS}")
  message(STATUS "netgen libs ${NETGEN_LIBRARY_DIR}")

  # Find zlib
  if(WIN32)
    add_definitions(-D_AFXDLL) # Needed by Netgen headers.
    find_package(ZLIB REQUIRED) # Find zlib headers since Netgen doesn't
    include_directories(${ZLIB_INCLUDE_DIRS})
    get_filename_component(ZLIB_LIB_PATH ${ZLIB_LIBRARIES} PATH)
  endif()
endif()

# Find ExodusII
if(ENABLE_EXODUS)
  find_library(EXO_LIB exoIIv2)
  find_file(EXO_HDR exodusII.h)
  get_filename_component(EXO_LIB_PATH ${EXO_LIB} PATH)
  get_filename_component(EXO_INCPATH ${EXO_HDR} PATH)
  list(APPEND LIBRARY_DIRS ${EXO_LIB_PATH})
  message(STATUS "exodus include path ${EXO_INCPATH}")
  message(STATUS "exodus libs ${EXO_LIB}")
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
  find_library(SIMMDL_LIB SimModel HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMADM_LIB SimAdvMeshing HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMDSC_LIB SimDiscrete HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMEXP_LIB SimExport HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMLCN_LIB SimLicense HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMMSH_LIB SimMeshing HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMMTL_LIB SimMeshTools HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMPMH_LIB SimPartitionedMesh HINTS ${SIMMETRIX_LIB_DIR})
  find_library(SIMPWR_LIB SimPartitionWrapper HINTS ${SIMMETRIX_LIB_DIR})
  set(SIMMETRIX_LIBS ${SIMMDL_LIB} ${SIMADM_LIB} ${SIMDSC_LIB} ${SIMEXP_LIB}
      ${SIMLCN_LIB} ${SIMMSH_LIB} ${SIMMTL_LIB} ${SIMPMH_LIB} ${SIMPWR_LIB})
  find_file(SIMMETRIX_HDR MeshSim.h HINTS ${SIMMETRIX_LIB_DIR}/../../include)
  get_filename_component(SIMMETRIX_INCPATH ${SIMMETRIX_HDR} PATH)
  get_filename_component(SIMMETRIX_LIBPATH ${SIMMDL_LIB} PATH)
  message(STATUS "Simmetrix include path ${SIMMETRIX_INCPATH}")
  message(STATUS "Simmetrix libs ${SIMMETRIX_LIBS}")
endif()

# Find boost
if(ENABLE_SIMMETRIX OR WIN32) # since WIN32 always needs boost
  if(WIN32)
    set(Boost_USE_STATIC_LIBS ON) # needed ON on WIN32. TODO: test on UNIX
  endif()
  set(Boost_USE_MULTITHREADED ON)
  set(Boost_USE_STATIC_RUNTIME OFF)
  find_package(Boost 1.45.0 REQUIRED COMPONENTS filesystem regex)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    list(APPEND LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
  endif()
endif()

## add include files to path
# TODO: Move this to targets.
include_directories(include)
#include_directories(AFTER ${MADLIB_INCPATH})
include_directories(AFTER ${GMSH_INCPATH})
include_directories(AFTER ${CGNS_INCPATH})
include_directories(AFTER ${METIS_INCPATH})
#include_directories(AFTER ${NETGEN_INCPATH})
include_directories(AFTER ${NETGEN_INCLUDE_DIRS})
include_directories(AFTER ${SIMMETRIX_INCPATH})
include_directories(AFTER ${HDF5_INCLUDE_DIRS})

# For some reason, older versions of cmake fail to find the local gmsh include
# files when compiling MAdLib header references.
# TODO: Is this still true in Gmsh 4.2+ and CMake 2.8+?
if(EXISTS ${GMSH_INCPATH})
  include_directories(AFTER ${GMSH_INCPATH}/../)
endif()

# Setting compile and link options #############################################

if(NOT WIN32) # These tests only work on GNU compilers.
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-std=c++14" COMPILER_SUPPORTS_CXX14)
  check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
  check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -DHAVE_TIMES")
    set(CMAKE_CXX_STANDARD 14)
  elseif(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DHAVE_TIMES")
    set(CMAKE_CXX_STANDARD 11)
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -DHAVE_TIMES")
  else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no
                    modern C++ support. Please use a different C++ compiler.")
  endif()
endif(NOT WIN32)

# Adding flags to Debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# Compiled objects output paths
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Turn off warnings  TODO: Remove as this gives appereance of clean code!
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
endif()

# Other flags for dependencies
if(ENABLE_SIMMETRIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_SIMMETRIX")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DHAVE_SIMMETRIX")
  message(STATUS "Using -DHAVE_SIMMETRIX")
endif()

# Setting library files ########################################################

# Setting include files
file(GLOB H_FILES include/*.H)
file(GLOB h_FILES include/*.h)
file(GLOB hpp_FILES include/*.hpp)
set(Nemosys_INC_FILES ${H_FILES} ${h_FILES} ${hpp_FILES})

# Setting r8 source
set(R8_SRCS
    src/math/r8lib.cpp
)

# Setting interp source
set(INTERP_SRCS
    src/math/baseInterp.C
    src/math/rbf_interp_nd.cpp
    src/math/rbfInterp.C
    src/spheres.C
)

# Setting Nemosys source
set(NEMOSYS_SRCS
    src/Drivers/NemDriver.C
    src/Drivers/ConversionDriver.C
    src/Drivers/InputGenDriver.C
    src/Drivers/MeshGenDriver.C
    src/Drivers/MeshQualityDriver.C
    src/Drivers/RefineDriver.C
    src/Drivers/RemeshDriver.C
    src/Drivers/RocPartCommGenDriver.C
#    src/Drivers/RocPrepDriver.C
#    src/Drivers/RocRestartDriver.C
    src/Drivers/TransferDriver.C

    src/InputGeneration/inputGen.C

    src/Integration/Cubature.C

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
    src/MeshPartitioning/meshStitcher.C

    src/MeshQuality/MeshQuality.C

    src/PatchRecovery/orthoPoly1D.C
    src/PatchRecovery/orthoPoly3D.C
    src/PatchRecovery/patchRecovery.C
    src/PatchRecovery/polyApprox.C

    src/Refinement/Refine.C

    src/SizeFieldGeneration/SizeFieldBase.C
    src/SizeFieldGeneration/GradSizeField.C
    src/SizeFieldGeneration/ValSizeField.C
    src/SizeFieldGeneration/Z2ErrorSizeField.C

    src/SolutionVerification/OrderOfAccuracy.C
#    src/SolutionVerification/RichardsonExtrapolation.C

    src/Transfer/TransferBase.C
    src/Transfer/FETransfer.C

    src/cgnsAnalyzer.C
    src/cgnsWriter.C
    src/gridTransfer.C
#    src/MAdLibInterface.C
#    src/meshPhys.C
#    src/netgenInterface.C
    src/rocstarCgns.C
    src/StlToVtk.C
    src/vtkAnalyzer.C
)

# Optional source files
if(ENABLE_EXODUS)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/Mesh/exoMesh.C
      src/MeshOperation/meshSrch.C
  )
endif()

if(ENABLE_EPIC)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/InputGeneration/ep16Prep.C
  )
endif()

if(ENABLE_SIMMETRIX)
  set(NEMOSYS_SRCS ${NEMOSYS_SRCS}
      src/MeshGeneration/simmetrixGen.C
      src/MeshGeneration/simmetrixParams.C
  )
endif()

set(UTIL_SRCS
    utils/Nemosys.C
    utils/cgns2msh.C
    utils/rocRemesh.C
    utils/rocSurfRemesh.C
    utils/xmlDump.C
    utils/rocStitchMesh.C
    utils/grid2gridTransfer.C
)
set(TUTR_SRCS
    utils/writeStructured.C
    tutorials/mshTransfer.C
    tutorials/dtkTutorial.C
)

# Building libraries ###########################################################

# Build r8 and rbf library
add_library(r8 ${R8_SRCS})

# Build interp library
add_library(interp ${INTERP_SRCS})
target_link_libraries(interp
    r8
    ANN
)

# Build Nemosys library
add_library(Nemosys ${NEMOSYS_SRCS})
target_link_libraries(Nemosys
    ${VTK_LIBRARIES}
    ${GMSH_LIB} TKernel
    ${CGNS_LIB}
    ${HDF5_LIBRARIES}
    ${Boost_LIBRARIES}
    interp
    MAdLib
)

# TODO: Appears unused.
#if(NOT ${BUILD_SHARED_LIBS})
#  add_definitions(-DSTATIC_LINK)
#endif()

if(ENABLE_MPI)
  target_link_libraries(Nemosys ${MPI_LIBRARIES})
endif()

if(ENABLE_EXODUS)
  target_link_libraries(Nemosys ${EXO_LIB})
  target_compile_definitions(Nemosys PRIVATE HAVE_EXODUSII)
endif()

if(ENABLE_EPIC)
  target_compile_definitions(Nemosys PRIVATE HAVE_EPIC)
endif()

if(ENABLE_SIMMETRIX)
  target_link_libraries(Nemosys ${SIMMETRIX_LIBS})
endif()

if(ENABLE_METIS)
  target_link_libraries(Nemosys ${METIS_LIB})
  target_compile_definitions(Nemosys PRIVATE HAVE_METIS)
endif()

if(ENABLE_NETGEN)
  target_link_libraries(Nemosys nglib)
  target_compile_definitions(Nemosys PRIVATE HAVE_NGEN)
endif()

# Install libraries ############################################################

# Headers
install(FILES ${Nemosys_INC_FILES}
    DESTINATION include)

# Libraries
install(TARGETS r8 interp Nemosys
    EXPORT Nemosys
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Export Nemosys as a CMake package
install(FILES NemosysConfig.cmake
    DESTINATION lib/cmake/Nemosys)
install(EXPORT Nemosys
    NAMESPACE Nemosys::
    DESTINATION lib/cmake/Nemosys)

# Setting utility files ########################################################

# Setting nemosysRun source
set(NEMOSYSRUN_SRCS utils/Nemosys.C)

if(ENABLE_BUILD_UTILS)
  # Setting xmlDump utility source
  set(XMLDUMP_SRCS utils/xmlDump.C)

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

  # Building DTK tutorial
  if(ENABLE_DTK)
    add_executable(dtkTutorial ${DTKTUTORIAL_SRC})
    set_target_properties(dtkTutorial PROPERTIES
        COMPILE_FLAGS ${dtk_headers}
        LINK_FLAGS ${dtk_libs})
    target_link_libraries(dtkTutorial Nemosys
        ${DTK_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    install(TARGETS dtkTutorial
        RUNTIME DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()
endif()

# Install utilities ############################################################

install(TARGETS nemosysRun
    RUNTIME DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

if(ENABLE_BUILD_UTILS)
  install(TARGETS xmlDump
                  cgns2msh
                  rocStitchMesh
                  rocRemesh
                  rocSurfRemesh
                  grid2gridTransfer
                  mshTransfer
      RUNTIME DESTINATION bin
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  if(ENABLE_DTK)
    install(TARGETS dtkTutorial
        RUNTIME DESTINATION bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif()
endif()

# Building Python wrapper ######################################################
if(ENABLE_PYTHON_BINDINGS)
  add_subdirectory(python)
endif()

# Building converters ##########################################################
if(ENABLE_EXODUS)
  add_subdirectory(converter)
endif()

# Building documentation #######################################################
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
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION share)
endif()

# Perform testing ##############################################################
if(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(testing)
endif(ENABLE_TESTING)

# Packaging ####################################################################
if(WIN32)
  list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/nemosysRun.exe")
  if(ENABLE_BUILD_UTILS)
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/xmlDump.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/cgns2msh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/rocStitchMesh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/rocRemesh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/rocSurfRemesh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/grid2gridTransfer.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/bin/mshTr.exe")
  endif()
  install(CODE "include(BundleUtilities)
  fixup_bundle(\"${CMAKE_EXE_TO_BUNDLE}\" \"\" \"${LIBRARY_DIRS}\")")
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION "bin")
  set(CMAKE_INSTALL_UCRT_LIBRARIES "TRUE")
  include(InstallRequiredSystemLibraries)
endif()

#set(CPACK_GENERATOR NSIS)
#set(CPACK_PACKAGE_NAME "MyLib")
#set(CPACK_PACKAGE_VENDOR "CMake.org")
#set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MyLib - CPack Component Installation Example")
#set(CPACK_PACKAGE_VERSION "1.0.0")
#set(CPACK_PACKAGE_VERSION_MAJOR "1")
#set(CPACK_PACKAGE_VERSION_MINOR "0")
#set(CPACK_PACKAGE_VERSION_PATCH "0")
#set(CPACK_PACKAGE_INSTALL_DIRECTORY "CPack Component Example")
#include(CPack)
