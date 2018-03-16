CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(Nemosys)
ENABLE_TESTING()

# this variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of Nemosys: e.g. "cmake
# -DDEFAULT=0 -DENABLE_DTK=1")
SET(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")

MACRO(opt OPTION HELP VALUE)
  OPTION(ENABLE_${OPTION} ${HELP} ${VALUE})
  SET(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
ENDMACRO(opt)

# writing remake scripts
FILE(WRITE ${CMAKE_BINARY_DIR}/makeAll
"#!/bin/bash
 rm -rf build
 mkdir build && cd build
 CMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}/madlib:${CMAKE_INSTALL_PREFIX}/gmsh:${CMAKE_INSTALL_PREFIX}/cgns:${CMAKE_INSTALL_PREFIX}/netgen cmake -DENABLE_BUILD_UTILS=ON -DENABLE_PYTHON_BINDINGS=ON -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} .. 
 make -j
 make install")
 FILE(WRITE ${CMAKE_BINARY_DIR}/makeNemOnly
"#!/bin/bash
 rm -rf build
 mkdir build && cd build
 CMAKE_PREFIX_PATH=${CMAKE_INSTALL_PREFIX}/madlib:${CMAKE_INSTALL_PREFIX}/gmsh:${CMAKE_INSTALL_PREFIX}/cgns:${CMAKE_INSTALL_PREFIX}/netgen cmake -DENABLE_BUILD_UTILS=ON -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} .. 
 make -j
 make install")
FILE(COPY ${CMAKE_BINARY_DIR}/makeAll DESTINATION ${CMAKE_SOURCE_DIR}
     FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
FILE(COPY ${CMAKE_BINARY_DIR}/makeNemOnly DESTINATION ${CMAKE_SOURCE_DIR}
     FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
     
# set to $DEFAULT to use by DEFAULT
opt(MPI "Enable MPI support" OFF)
opt(DTK "Enable DTK extensions" OFF)
opt(EXODUS "Enable EXODUS II extensions" OFF)
opt(PYTHON_BINDINGS "Enable Python bindings" OFF)
opt(TESTING "Enable testing" ON)
opt(BUILD_UTILS "Build utilities" OFF)

# set up version
SET(NEMOSYS_MAJOR_VERSION 1)
SET(NEMOSYS_MINOR_VERSION 2)
SET(NEMOSYS_PATCH_VERSION 0)
SET(NEMOSYS_EXTRA_VERSION "" CACHE STRING "NEMoSys extra version string")
SET(NEMOSYS_VERSION "${NEMOSYS_MAJOR_VERSION}.${NEMOSYS_MINOR_VERSION}")
SET(NEMOSYS_VERSION "${NEMOSYS_VERSION}.${NEMOSYS_PATCH_VERSION}${NEMOSYS_EXTRA_VERSION}")
SET(NEMOSYS_SHORT_LICENSE "GNU Lesser General Public License (LGPL)")

# finding libraries ###############################################################
IF (ENABLE_MPI)
  FIND_PACKAGE(MPI)
  IF(MPI_CXX_COMPILER)
    SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  ENDIF()
ENDIF (ENABLE_MPI)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_LIBRARY(MADLIB_LIB MAdLib HINTS ${CMAKE_INSTALL_PREFIX}/madlib/lib)
FIND_FILE(MADLIB_HDR MAdLib.h HINTS ${CMAKE_INSTALL_PREFIX}/madlib/include/MAdLib)
GET_FILENAME_COMPONENT(MADLIB_INCPATH ${MADLIB_HDR} PATH)
MESSAGE(STATUS "madlib include path: ${MADLIB_INCPATH}")

FIND_LIBRARY(GMSH_LIB Gmsh HINTS ${CMAKE_INSTALL_PREFIX}/gmsh/lib)
FIND_FILE(GMSH_HDR Gmsh.h HINTS ${CMAKE_INSTALL_PREFIX}/gmsh/include/gmsh)
GET_FILENAME_COMPONENT(GMSH_INCPATH ${GMSH_HDR} PATH)
MESSAGE(STATUS "gmsh include path ${GMSH_INCPATH}")

FIND_LIBRARY(CGNS_LIB cgns HINTS ${CMAKE_INSTALL_PREFIX}/cgns/lib)
FIND_FILE(CGNS_HDR cgnslib.h HINTS ${CMAKE_INSTALL_PREFIX}/include)
GET_FILENAME_COMPONENT(CGNS_INCPATH ${CGNS_HDR} PATH)
MESSAGE(STATUS "cgns include path ${CGNS_INCPATH}")

FIND_LIBRARY(METIS_LIB metis)
FIND_FILE(METIS_HDR metis.h)
GET_FILENAME_COMPONENT(METIS_INCPATH ${METIS_HDR} PATH)

FIND_LIBRARY(NETGEN_LIB nglib HINTS ${CMAKE_INSTALL_PREFIX}/netgen/lib/)
FIND_FILE(NETGEN_HDR nglib.h HINTS ${CMAKE_INSTALL_PREFIX}/netgen/include/include)
GET_FILENAME_COMPONENT(NETGEN_INCPATH ${NETGEN_HDR} PATH)
MESSAGE(STATUS "netgen include path ${NETGEN_INCPATH}")


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

## add include files to path
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(AFTER ${MADLIB_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${GMSH_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${CGNS_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${METIS_INCPATH})
INCLUDE_DIRECTORIES(AFTER ${NETGEN_INCPATH})

# For some reason, older versions of cmake fail to find the local gmsh include
# files when compiling madlib header references.
IF(EXISTS ${GMSH_INCPATH})
  INCLUDE_DIRECTORIES(AFTER ${GMSH_INCPATH}/../)
ENDIF()

# Setting compile and link options ############################################
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

### adding flags to debug
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

# Compiled objects output paths
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/netgen/lib:${CMAKE_INSTALL_PREFIX}/gmsh/lib:${CMAKE_INSTALL_PREFIX}/madlib/lib:${CMAKE_INSTALL_PREFIX}/Nemosys/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# turning off warnings
IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
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
                 src/netgenInterface.C src/MeshGeneration/meshGen.C
                 src/Transfer/TransferBase.C  src/Transfer/FETransfer.C
                 src/SizeFieldGeneration/SizeFieldBase.C
                 src/SizeFieldGeneration/GradSizeField.C
                 src/SizeFieldGeneration/ValSizeField.C
                 src/SizeFieldGeneration/Z2ErrorSizeField.C
                 src/Refinement/Refine.C src/MeshQuality/MeshQuality.C
                 src/Drivers/NemDriver.C src/Drivers/TransferDriver.C
                 src/Drivers/RefineDriver.C src/Drivers/MeshGenDriver.C
                 src/Drivers/MeshQualityDriver.C src/StlToVtk.C src/Integration/Cubature.C
                 src/PatchRecovery/orthoPoly1D.C src/PatchRecovery/orthoPoly3D.C
                 src/PatchRecovery/polyApprox.C src/PatchRecovery/patchRecovery.C)
SET(ALTNEMOSYS_SRCS src/vtkAnalyzer.C src/cgnsAnalyzer.C src/meshPartitioner.C
                src/cgnsWriter.C src/gridTransfer.C)
SET(ROCSTAR_SRCS src/rocstarCgns.C)
SET(UTIL_SRCS utils/Nemosys.C utils/cgns2msh.C utils/rocRemesh utils/rocSurfRemesh.C)
SET(TUTR_SRCS utils/writeStructured.C)
SET(R8_SRCS src/math/r8lib.cpp)
SET(SPHERE_SRCS src/spheres.C)
SET(INTERP_SRCS src/math/rbf_interp_nd.cpp src/math/rbfInterp.C src/math/baseInterp.C)
SET(ALL_SRCS ${NEMOSYS_SRCS} ${ALTNEMOSYS_SRCS} ${UTIL_SRCS} ${TUTR_SRCS}
             ${R8_SRCS} ${INTERP_SRCS} ${ROCSTAR_SRCS} ${SPHERE_SRCS})
SET_SOURCE_FILES_PROPERTIES(${ALL_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC -O3 ${MPI_CXX_COMPILE_FLAGS}" )
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
  ADD_LIBRARY(interp SHARED ${INTERP_SRCS} ${SPHERE_SRCS})
ELSE()
  ADD_LIBRARY(interp ${R8_SRCS})
ENDIF()
TARGET_LINK_LIBRARIES(interp r8)

# build nemosys lib
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(Nemosys SHARED ${NEMOSYS_SRCS})
ELSE()
  ADD_LIBRARY(Nemosys ${NEMOSYS_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()
TARGET_LINK_LIBRARIES(Nemosys ${MADLIB_LIB} ${NETGEN_LIB} ${GMSH_LIB} ${VTK_LIBRARIES})

# build old nemosys lib
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(AltNemosys SHARED ${ALTNEMOSYS_SRCS})
ELSE()
  ADD_LIBRARY(AltNemosys ${ALTNEMOSYS_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()
TARGET_LINK_LIBRARIES(AltNemosys ${MADLIB_LIB} ${GMSH_LIB} ${VTK_LIBRARIES} ${CGNS_LIB} ${METIS_LIB} interp)

# Rocstar support library
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(Rocs SHARED ${ROCSTAR_SRCS})
ELSE()
  ADD_LIBRARY(Rocs ${ROCSTAR_SRCS})
ENDIF()

INSTALL(TARGETS r8 interp Nemosys AltNemosys Rocs
        LIBRARY DESTINATION Nemosys/lib PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

# Builiding utilities #########################################################

SET(CMAKE_CXX_LINK_EXECUTABLE_OLD ${CMAKE_CXX_LINK_EXECUTABLE})
SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lfltk_images -lfltk -lfltk_gl -lfltk_forms -lSM -lICE -lX11 -lXext -lm -ljpeg -lz -lpng -lGLU -lGL -llapack -lgmp -lblas -lcgns -lhdf5_cpp")

# Building Nemosys Core utility ###########################################
ADD_EXECUTABLE(nemosysRun utils/Nemosys.C)
IF(VTK_LIBRARIES)
  TARGET_LINK_LIBRARIES(nemosysRun Nemosys)
ELSE()
  TARGET_LINK_LIBRARIES(nemosysRun Nemosys vtkHybrid vtkWidgets)
ENDIF()
INSTALL(TARGETS nemosysRun RUNTIME DESTINATION Nemosys/bin
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
  TARGET_LINK_LIBRARIES(cgns2msh AltNemosys)

  # Rocstart volume remesh utility
  ADD_EXECUTABLE(rocRemesh utils/rocRemesh.C)
  TARGET_LINK_LIBRARIES(rocRemesh AltNemosys)

  # Rocstar surface remesh utility
  ADD_EXECUTABLE(rocSurfRemesh utils/rocSurfRemesh.C)
  TARGET_LINK_LIBRARIES(rocSurfRemesh AltNemosys Rocs)

  # Rocstar surface remesh utility
  ADD_EXECUTABLE(grid2gridTransfer utils/grid2gridTransfer.C)
  TARGET_LINK_LIBRARIES(grid2gridTransfer AltNemosys Rocs)

  # Building tutorials ###########################################################
  SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE_OLD} -lhdf5_cpp")

  # meshTransfer tutorial
  ADD_EXECUTABLE(mshTransfer tutorials/mshTransfer.C)
  TARGET_LINK_LIBRARIES(mshTransfer AltNemosys)
  
  INSTALL(TARGETS xmlDump cgns2msh rocRemesh rocSurfRemesh grid2gridTransfer mshTransfer 
          RUNTIME DESTINATION Nemosys/bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
          GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  # DTK tutorial
  IF (ENABLE_DTK)
    ADD_EXECUTABLE(dtkTutorial tutorials/dtkTutorial.C)
    SET_TARGET_PROPERTIES(dtkTutorial PROPERTIES COMPILE_FLAGS ${dtk_headers} LINK_FLAGS ${dtk_libs})
    TARGET_LINK_LIBRARIES(dtkTutorial AltNemosys Rocs
        ${DTK_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    INSTALL(TARGETS dtkTutorial RUNTIME DESTINATION Nemosys/bin
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
SET(CMAKE_CXX_LINK_EXECUTABLE ${CMAKE_CXX_LINK_EXECUTABLE_OLD})
IF (ENABLE_TESTING)
    # adding gtest subdir with its own cmake file
    ADD_SUBDIRECTORY(testing)
    ADD_SUBDIRECTORY(testing/lib/gtest-1.6.0)
    # include gtest library. gtest_SOURCE_DIR is availble from subdir addition
    INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}/include/gtest ${gtest_SOURCE_DIR})
    SET(OLD_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing/test_scripts)
    ADD_EXECUTABLE(runCubatureInterpTest testing/test_scripts/testCubature.C)
    ADD_EXECUTABLE(runConversionTest testing/test_scripts/testConversion.C)
    ADD_EXECUTABLE(runOrthoPolyTest testing/test_scripts/testOrthoPoly.C)
    ADD_EXECUTABLE(runPatchRecoveryTest testing/test_scripts/testPatchRecovery.C)
    ADD_EXECUTABLE(runTransferTest testing/test_scripts/testTransfer.C)
    ADD_EXECUTABLE(runMeshGenTest testing/test_scripts/testMeshGen.C)
    TARGET_LINK_LIBRARIES(runCubatureInterpTest gtest gtest_main Nemosys)
    TARGET_LINK_LIBRARIES(runConversionTest gtest gtest_main Nemosys)
    TARGET_LINK_LIBRARIES(runOrthoPolyTest gtest gtest_main Nemosys)
    TARGET_LINK_LIBRARIES(runPatchRecoveryTest gtest gtest_main Nemosys)
    TARGET_LINK_LIBRARIES(runTransferTest gtest gtest_main Nemosys)
    TARGET_LINK_LIBRARIES(runMeshGenTest gtest gtest_main Nemosys)
    SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_RUNTIME_OUTPUT_DIRECTORY})
ENDIF(ENABLE_TESTING)
