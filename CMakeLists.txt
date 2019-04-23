CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(Nemosys)
ENABLE_TESTING()

# this variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of Nemosys: e.g. "cmake
# -DDEFAULT=0 -DENABLE_DTK=1")
SET(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")
# creating code compilation database
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )

MACRO(opt OPTION HELP VALUE)
  OPTION(ENABLE_${OPTION} ${HELP} ${VALUE})
  SET(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
ENDMACRO(opt)

## writing remake scripts
if(NOT WIN32)
  FILE(WRITE ${CMAKE_BINARY_DIR}/makeAll
  "#!/bin/bash
   rm -rf build
   mkdir build && cd build
   CMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}/madlib:${CMAKE_INSTALL_PREFIX}/gmsh:${CMAKE_INSTALL_PREFIX}/netgen cmake -DENABLE_BUILD_UTILS=ON -DENABLE_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ..
   make -j
   make install")
  FILE(WRITE ${CMAKE_BINARY_DIR}/makeNemOnly
  "#!/bin/bash
   rm -rf build
   mkdir build && cd build
   CMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}/madlib:${CMAKE_INSTALL_PREFIX}/gmsh:${CMAKE_INSTALL_PREFIX}/netgen cmake -DENABLE_BUILD_UTILS=ON -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ..
   make -j
   make install")
  FILE(COPY ${CMAKE_BINARY_DIR}/makeAll DESTINATION ${CMAKE_SOURCE_DIR}
       FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
  FILE(COPY ${CMAKE_BINARY_DIR}/makeNemOnly DESTINATION ${CMAKE_SOURCE_DIR}
       FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
endif()

# set to $DEFAULT to use by DEFAULT
opt(MPI "Enable MPI support" OFF)
opt(DTK "Enable DTK extensions" OFF)
opt(EXODUS "Enable EXODUS II extensions" OFF)
opt(EPIC "Enable EPIC preprocessor" OFF)
opt(PYTHON_BINDINGS "Enable Python bindings" OFF)
opt(TESTING "Enable testing" ON)
opt(BUILD_UTILS "Build utilities" OFF)
opt(METIS "Enable Metis interface" ON)
opt(NETGEN "Enable Netgen interface" ON)
opt(SYMMX "Enable SYMMETRIX Meshing engine" OFF)

# set up version
SET(NEMOSYS_MAJOR_VERSION 0)
SET(NEMOSYS_MINOR_VERSION 28)
SET(NEMOSYS_PATCH_VERSION 0)
SET(NEMOSYS_EXTRA_VERSION "" CACHE STRING "NEMoSys extra version string")
SET(NEMOSYS_VERSION "${NEMOSYS_MAJOR_VERSION}.${NEMOSYS_MINOR_VERSION}")
SET(NEMOSYS_VERSION "${NEMOSYS_VERSION}.${NEMOSYS_PATCH_VERSION}${NEMOSYS_EXTRA_VERSION}")
SET(NEMOSYS_SHORT_LICENSE "GNU Lesser General Public License (LGPL)")

# finding libraries ###############################################################
IF (ENABLE_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  # include mpi headers
  include_directories(${MPI_INCLUDE_PATH})
  IF (MPI_LIBRARY)
    SET(MPI_LIBRARIES ${MPI_LIBRARY})
  ENDIF (MPI_LIBRARY)
  IF (MPI_EXTRA_LIBRARY)
    SET(MPI_LIBRARIES ${MPI_LIBRARIES} "${MPI_EXTRA_LIBRARY}")
  ENDIF (MPI_EXTRA_LIBRARY)
  # Needed for mpich 2
  ADD_DEFINITIONS("-DMPICH_IGNORE_CXX_SEEK")
ENDIF (ENABLE_MPI)

if(WIN32)
  SET(VTK_BUILD_SHARED_LIBS OFF)
endif()
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_UPPER)
get_target_property(VTK_LIB_DIR vtksys IMPORTED_LOCATION_${CMAKE_BUILD_TYPE_UPPER})
GET_FILENAME_COMPONENT(VTK_LIB_DIR ${VTK_LIB_DIR} PATH)
LIST(APPEND LIBRARY_DIRS ${VTK_LIB_DIR})
MESSAGE(STATUS "vtk include path: ${VTK_INCLUDE_DIRS}")
message(STATUS "vtk libs ${VTK_LIB_DIR}")

add_subdirectory(contrib/madlib-1.3.0)

FIND_LIBRARY(GMSH_LIB gmsh HINTS ${CMAKE_INSTALL_PREFIX}/gmsh/lib)
FIND_FILE(GMSH_HDR gmsh.h HINTS ${CMAKE_INSTALL_PREFIX}/gmsh/include)
GET_FILENAME_COMPONENT(GMSH_INCPATH ${GMSH_HDR} PATH)
MESSAGE(STATUS "gmsh include path ${GMSH_INCPATH}")

if(WIN32)
  SET(HDF5_USE_STATIC_LIBRARIES ON) # Needed on WIN32. TODO: Test on UNIX
endif()
FIND_PACKAGE(HDF5 REQUIRED)
GET_FILENAME_COMPONENT(HDF5_LIB_PATH ${HDF5_LIBRARIES} PATH)
LIST(APPEND LIBRARY_DIRS ${HDF5_LIB_PATH})
message(STATUS "hdf5 include path ${HDF5_INCLUDE_DIRS}")
message(STATUS "hdf5 libs ${HDF5_LIBRARIES}")

FIND_LIBRARY(CGNS_LIB cgns HINTS ${CMAKE_INSTALL_PREFIX}/cgns/lib)
FIND_FILE(CGNS_HDR cgnslib.h HINTS ${CMAKE_INSTALL_PREFIX}/cgns/include)
GET_FILENAME_COMPONENT(CGNS_INCPATH ${CGNS_HDR} PATH)
GET_FILENAME_COMPONENT(CGNS_LIB_PATH ${CGNS_LIB} PATH)
LIST(APPEND LIBRARY_DIRS ${CGNS_LIB_PATH})
MESSAGE(STATUS "cgns include path ${CGNS_INCPATH}")
MESSAGE(STATUS "cgns libs ${CGNS_LIB}")

IF (ENABLE_METIS)
  FIND_LIBRARY(METIS_LIB metis HINTS ${CMAKE_INSTALL_PREFIX}/metis/lib)
  FIND_FILE(METIS_HDR metis.h HINTS ${CMAKE_INSTALL_PREFIX}/metis/include)
  GET_FILENAME_COMPONENT(METIS_INCPATH ${METIS_HDR} PATH)
  GET_FILENAME_COMPONENT(METIS_LIB_PATH ${METIS_LIB} PATH)
  LIST(APPEND LIBRARY_DIRS ${METIS_LIB_PATH})
  MESSAGE(STATUS "metis include path ${METIS_INCPATH}")
  MESSAGE(STATUS "metis libs ${METIS_LIB}")
  ADD_DEFINITIONS(-DHAVE_METIS)
ENDIF (ENABLE_METIS)

IF (ENABLE_NETGEN)
  FIND_LIBRARY(NETGEN_LIB nglib HINTS ${CMAKE_INSTALL_PREFIX}/netgen/lib)
  FIND_FILE(NETGEN_HDR nglib.h HINTS ${CMAKE_INSTALL_PREFIX}/netgen/include)
  GET_FILENAME_COMPONENT(NETGEN_INCPATH ${NETGEN_HDR} PATH)
  GET_FILENAME_COMPONENT(NETGEN_LIB_PATH ${NETGEN_LIB} PATH)
  LIST(APPEND LIBRARY_DIRS ${NETGEN_LIB_PATH})
  MESSAGE(STATUS "netgen include path ${NETGEN_INCPATH}")
  MESSAGE(STATUS "netgen libs ${NETGEN_LIB}")
  ADD_DEFINITIONS(-DHAVE_NGEN)
  if(WIN32)
    ADD_DEFINITIONS(-D_AFXDLL) # Needed by Netgen headers.
  	FIND_PACKAGE(ZLIB REQUIRED) # Find zlib headers since Netgen doesn't
  	include_directories(${ZLIB_INCLUDE_DIRS})
    GET_FILENAME_COMPONENT(ZLIB_LIB_PATH ${ZLIB_LIBRARIES} PATH)
  endif(WIN32)
ENDIF (ENABLE_NETGEN)

IF (ENABLE_DTK)
  FIND_PACKAGE(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})
  IF(NOT Trilinos_FOUND)
    MESSAGE(FATAL_ERROR "Could not find Trilinos!")
  ENDIF()
  MESSAGE (Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS})
  LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS})
  LINK_DIRECTORIES(${Trilinos_TPL_LIBRARY_DIRS})
ENDIF()

# import dtk related variables
IF (ENABLE_DTK)
  INCLUDE(dtk_stuff.txt)
ENDIF ()

# symmetrix
IF (ENABLE_SYMMX)
  FIND_LIBRARY(SYMMDL_LIB SimModel HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMADM_LIB SimAdvMeshing HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMDSC_LIB SimDiscrete HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMEXP_LIB SimExport HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMLIC_LIB SimLicense HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMMSH_LIB SimMeshing HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMMTL_LIB SimMeshTools HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMPMS_LIB SimPartitionedMesh HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  FIND_LIBRARY(SYMPRP_LIB SimPartitionWrapper HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/lib/x64_rhel7_gcc48)
  SET(SYMMX_LIBS ${SYMADM_LIB} ${SYMEXP_LIB} ${SYMDSC_LIB} ${SYMPMS_LIB} ${SYMPRP_LIB} ${SYMMSH_LIB} ${SYMMTL_LIB} ${SYMMDL_LIB} ${SYMLIC_LIB})
  FIND_FILE(SYMMX_HDR MeshSim.h HINTS /home/users/snatesh/snatesh-storage/SYMMETRIX/11.0-180502/include)
  GET_FILENAME_COMPONENT(SYMMX_INCPATH ${SYMMX_HDR} PATH)
  GET_FILENAME_COMPONENT(SYMMX_LIBPATH ${SYMMDL_LIB} PATH)
  MESSAGE(STATUS "SYMMETRIX library path ${SYMMX_LIBPATH}")
  MESSAGE(STATUS "SYMMETRIX libraries ${SYMMX_LIBS}")
  MESSAGE(STATUS "SYMMETRIX include path ${SYMMX_INCPATH}")
ENDIF()

IF(ENABLE_SYMMX OR WIN32) # since WIN32 always needs boost
  # boost
  SET(Boost_USE_STATIC_LIBS ON) # needed ON on WIN32. TODO: test on UNIX
  SET(Boost_USE_MULTITHREADED ON)
  SET(Boost_USE_STATIC_RUNTIME OFF)
  FIND_PACKAGE(Boost 1.45.0 REQUIRED COMPONENTS filesystem regex)
  IF(Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
    LIST(APPEND LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
  ENDIF()
ENDIF()

## add include files to path
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(AFTER ${MADLIB_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${GMSH_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${CGNS_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${METIS_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${NETGEN_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${SYMMX_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${HDF5_INCLUDE_DIRS})

# For some reason, older versions of cmake fail to find the local gmsh include
# files when compiling madlib header references.
IF(EXISTS ${GMSH_INCPATH})
  INCLUDE_DIRECTORIES(AFTER ${GMSH_INCPATH}/../)
ENDIF()

# Setting compile and link options ############################################
IF(NOT WIN32) # These tests only work on GNU compilers.
  INCLUDE(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  IF(COMPILER_SUPPORTS_CXX14)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -DHAVE_TIMES")
  ELSEIF(COMPILER_SUPPORTS_CXX11)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DHAVE_TIMES")
  ELSEIF(COMPILER_SUPPORTS_CXX0X)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -DHAVE_TIMES")
  ELSE()
    MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has
                    no modern C++ support. Please use a different C++ compiler.")
  ENDIF()
ENDIF(NOT WIN32)

# adding flags to debug
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# compiled objects output paths
#SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# setting RPATH
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/netgen/lib:${CMAKE_INSTALL_PREFIX}/gmsh/lib:${CMAKE_INSTALL_PREFIX}/madlib/lib:${CMAKE_INSTALL_PREFIX}/Nemosys/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# turning off warnings
IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
ENDIF()

# additional compiler flags
IF (ENABLE_EXODUS)
  FIND_LIBRARY(EXO_LIB exodus)
  FIND_FILE(EXO_HDR exodusII.h)
  GET_FILENAME_COMPONENT(EXO_INCPATH ${EXO_HDR} PATH)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_EXODUSII")
  SET(NEMOSYS_LNK ${EXO_LIB})
ENDIF()

IF (ENABLE_EPIC)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_EPIC")
ENDIF()


# other flags for dependencies
IF (ENABLE_SYMMX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_SYMMX")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DHAVE_SYMMX")
  MESSAGE(STATUS "Using -DHAVE_SYMMX")
ENDIF()

# Setting files   ###############################################################
# Install include files
FILE(GLOB H_FILES include/*.H)
FILE(GLOB h_FILES include/*.h)
FILE(GLOB hpp_FILES include/*.hpp)
SET(Nemosys_INC_FILES ${H_FILES} ${h_FILES} ${hpp_FILES})
INSTALL(FILES ${Nemosys_INC_FILES} DESTINATION Nemosys/include)

# Setting compile and link flags
SET(NEMOSYS_SRCS src/Mesh/meshBase.C src/Mesh/vtkMesh.C
                 src/MeshGeneration/meshGen.C
                 src/MeshGeneration/netgenGen.C src/MeshGeneration/netgenParams.C
                 src/Transfer/TransferBase.C  src/Transfer/FETransfer.C
                 src/SizeFieldGeneration/SizeFieldBase.C
                 src/SizeFieldGeneration/GradSizeField.C
                 src/SizeFieldGeneration/ValSizeField.C
                 src/SizeFieldGeneration/Z2ErrorSizeField.C
                 src/Refinement/Refine.C src/MeshQuality/MeshQuality.C
                 src/Drivers/NemDriver.C src/Drivers/TransferDriver.C
                 src/Drivers/RefineDriver.C src/Drivers/RemeshDriver.C
                 src/Drivers/MeshGenDriver.C src/Drivers/MeshQualityDriver.C
                 src/StlToVtk.C src/Integration/Cubature.C
                 src/Drivers/RocPartCommGenDriver.C
                 src/PatchRecovery/orthoPoly1D.C src/PatchRecovery/orthoPoly3D.C
                 src/PatchRecovery/polyApprox.C src/PatchRecovery/patchRecovery.C
                 src/Mesh/pntMesh.C src/Drivers/ConversionDriver.C
                 src/SolutionVerification/OrderOfAccuracy.C src/vtkAnalyzer.C
                 src/cgnsAnalyzer.C src/rocstarCgns.C
                 src/MeshPartitioning/meshPartitioner.C
                 src/MeshPartitioning/meshStitcher.C
                 src/cgnsWriter.C src/gridTransfer.C
                 src/Mesh/cobalt.C src/Mesh/patran.C
                 src/MeshOperation/meshSrch.C
                 src/InputGeneration/inputGen.C
                 src/Drivers/InputGenDriver.C)

# optional source files
IF (ENABLE_EXODUS)
  SET(NEMOSYS_SRCS ${NEMOSYS_SRCS} src/Mesh/exoMesh.C )
ENDIF()

IF (ENABLE_EPIC)
  SET(NEMOSYS_SRCS ${NEMOSYS_SRCS} src/InputGeneration/ep16Prep.C )
ENDIF()

SET(ALTNEMOSYS_SRCS src/vtkAnalyzer.C src/cgnsAnalyzer.C src/meshPartitioner.C
                    src/cgnsWriter.C src/gridTransfer.C)
SET(ROCSTAR_SRCS src/rocstarCgns.C)
SET(UTIL_SRCS utils/Nemosys.C utils/cgns2msh.C utils/rocRemesh.C utils/rocSurfRemesh.C
              utils/xmlDump.C utils/rocStitchMesh.C utils/grid2gridTransfer.C)

SET(TUTR_SRCS utils/writeStructured.C tutorials/mshTransfer.C tutorials/dtkTutorial.C)
SET(R8_SRCS src/math/r8lib.cpp)
SET(INTERP_SRCS src/math/rbf_interp_nd.cpp src/math/rbfInterp.C src/math/baseInterp.C src/spheres.C)
IF(ENABLE_SYMMX)
  SET(SYMMX_SRCS src/MeshGeneration/symmxGen.C src/MeshGeneration/symmxParams.C)
  SET(NEMOSYS_SRCS ${NEMOSYS_SRCS} ${SYMMX_SRCS})
ENDIF()
SET(ALL_SRCS ${NEMOSYS_SRCS} ${ALTNEMOSYS_SRCS} ${UTIL_SRCS} ${TUTR_SRCS}
             ${R8_SRCS} ${INTERP_SRCS} ${SYMMX_SRCS})
#SET_SOURCE_FILES_PROPERTIES(${ALL_SRCS} APPEND PROPERTIES COMPILE_FLAGS "-fPIC -O3 ${MPI_CXX_COMPILE_FLAGS}" )
ADD_DEFINITIONS(-D_HAVE_GMSH_)


# Builiding libraries #########################################################

# Build r8 and rbf library
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(r8 SHARED ${R8_SRCS})
ELSE()
  ADD_LIBRARY(r8 ${R8_SRCS})
ENDIF()

# build interp lib
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(interp SHARED ${INTERP_SRCS})
ELSE()
  ADD_LIBRARY(interp ${R8_SRCS})
ENDIF()
TARGET_LINK_LIBRARIES(interp r8 ANN)

# build nemosys lib
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(Nemosys SHARED ${NEMOSYS_SRCS})
ELSE()
  ADD_LIBRARY(Nemosys ${NEMOSYS_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()
TARGET_LINK_LIBRARIES(Nemosys ${VTK_LIBRARIES} ${NETGEN_LIB} ${GMSH_LIB} ${SYMMX_LIBS}
                              ${CGNS_LIB} ${HDF5_LIBRARIES} ${METIS_LIB} ${Boost_LIBRARIES}
                              ${NEMOSYS_LNK}
                              interp MAdLib)

INSTALL(TARGETS r8 interp Nemosys
        RUNTIME DESTINATION Nemosys/bin
        LIBRARY DESTINATION Nemosys/lib
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)  

# Builiding utilities #########################################################

#SET(CMAKE_CXX_LINK_EXECUTABLE_OLD ${CMAKE_CXX_LINK_EXECUTABLE})
#SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lfltk_images -lfltk -lfltk_gl -lfltk_forms -lSM -lICE -lX11 -lXext -lm -ljpeg -lz -lpng -lGLU -lGL -llapack -lgmp -lblas -lcgns -lhdf5_cpp")

# Building Nemosys Core utility ###########################################
ADD_EXECUTABLE(nemosysRun utils/Nemosys.C)
IF(ENABLE_MPI)
  TARGET_LINK_LIBRARIES(nemosysRun Nemosys ${MPI_LIBRARIES})
ELSE()
  TARGET_LINK_LIBRARIES(nemosysRun Nemosys)
ENDIF()

INSTALL(TARGETS nemosysRun
        RUNTIME DESTINATION Nemosys/bin
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

IF(ENABLE_BUILD_UTILS)

  # xmlDump utility
  ADD_EXECUTABLE(xmlDump MACOSX_BUNDLE utils/xmlDump.C)
  if(VTK_LIBRARIES)
    TARGET_LINK_LIBRARIES(xmlDump Nemosys)
  else()
    TARGET_LINK_LIBRARIES(xmlDump Nemosys vtkHybrid vtkWidgets)
  ENDIF()

  # cgns2gmsh utility
  ADD_EXECUTABLE(cgns2msh utils/cgns2msh.C)
  TARGET_LINK_LIBRARIES(cgns2msh Nemosys)

  # stitch cgns utility
  ADD_EXECUTABLE(rocStitchMesh utils/rocStitchMesh.C)
  TARGET_LINK_LIBRARIES(rocStitchMesh Nemosys)

  # Rocstart volume remesh utility
  ADD_EXECUTABLE(rocRemesh utils/rocRemesh.C)
  TARGET_LINK_LIBRARIES(rocRemesh Nemosys)

  ADD_EXECUTABLE(rocSurfRemesh utils/rocSurfRemesh.C)
  TARGET_LINK_LIBRARIES(rocSurfRemesh Nemosys)

  # Rocstar surface remesh utility
  ADD_EXECUTABLE(grid2gridTransfer utils/grid2gridTransfer.C)
  TARGET_LINK_LIBRARIES(grid2gridTransfer Nemosys)

  # Building tutorials ###########################################################
  #SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE_OLD} -lhdf5_cpp")

  # meshTransfer tutorial
  ADD_EXECUTABLE(mshTransfer tutorials/mshTransfer.C)
  TARGET_LINK_LIBRARIES(mshTransfer Nemosys)

  INSTALL(TARGETS xmlDump cgns2msh rocStitchMesh rocRemesh rocSurfRemesh grid2gridTransfer mshTransfer
          RUNTIME DESTINATION Nemosys/bin
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
          GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  # DTK tutorial
  IF (ENABLE_DTK)
    ADD_EXECUTABLE(dtkTutorial tutorials/dtkTutorial.C)
    SET_TARGET_PROPERTIES(dtkTutorial PROPERTIES COMPILE_FLAGS ${dtk_headers} LINK_FLAGS ${dtk_libs})
    TARGET_LINK_LIBRARIES(dtkTutorial Nemosys
        ${DTK_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    INSTALL(TARGETS dtkTutorial
            RUNTIME DESTINATION Nemosys/bin
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  ENDIF()

ENDIF(ENABLE_BUILD_UTILS)

# Building python wrapper #####################################################
IF (ENABLE_PYTHON_BINDINGS)
  ADD_SUBDIRECTORY(python)
  file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/Nemosys/python/lib)
ENDIF()

# Building converters ##########################################################
IF (ENABLE_EXODUS)
  ADD_SUBDIRECTORY(converter)
ENDIF()

# Builiding documentation #####################################################
OPTION(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" OFF)
IF(BUILD_DOCUMENTATION)
  FIND_PACKAGE(Doxygen)
  IF(NOT DOXYGEN_FOUND)
    MESSAGE(FATAL_ERROR "Doxygen is needed to build the documentation.")
  ENDIF()
  SET( doxyfile_in          ${CMAKE_CURRENT_SOURCE_DIR}/doc/Nemosys.cfg.in )
  SET( doxyfile             ${CMAKE_CURRENT_BINARY_DIR}/Nemosys.cfg        )
  SET( documentation        ${CMAKE_CURRENT_BINARY_DIR}/docs                 )
  CONFIGURE_FILE( ${doxyfile_in} ${doxyfile} @ONLY )
  ADD_CUSTOM_COMMAND( OUTPUT ${documentation}
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
                      MAIN_DEPENDENCY ${doxyfile} ${doxyfile_in}
                      DEPENDS Nemosys ${doxy_extra_files}
                      COMMENT "Generating documentation")
  ADD_CUSTOM_TARGET( docs DEPENDS ${documentation} )
  INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION Nemosys/share )
ENDIF()

# Perform testing ############################################################
#SET(CMAKE_CXX_LINK_EXECUTABLE ${CMAKE_CXX_LINK_EXECUTABLE_OLD})
IF (ENABLE_TESTING)
  ADD_SUBDIRECTORY(testing)
ENDIF(ENABLE_TESTING)

if(WIN32)
  list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/nemosysRun.exe")
  if(ENABLE_BUILD_UTILS)
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/xmlDump.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/cgns2msh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/rocStitchMesh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/rocRemesh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/rocSurfRemesh.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/grid2gridTransfer.exe")
    list(APPEND CMAKE_EXE_TO_BUNDLE "${CMAKE_INSTALL_PREFIX}/Nemosys/bin/mshTr.exe")
  endif()
  install(CODE "include(BundleUtilities)
  fixup_bundle(\"${CMAKE_EXE_TO_BUNDLE}\" \"\" \"${LIBRARY_DIRS}\")")
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION "Nemosys/bin")
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
