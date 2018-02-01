CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(Nemosys)
ENABLE_TESTING()
# this variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of Nemosys: e.g. "cmake
# -DDEFAULT=0 -DENABLE_DTK=1")
set(DEFAULT ON CACHE INTERNAL "Default value for enabled-by-default options")

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro(opt)

# set to $DEFAULT to use by DEFAULT
opt(MPI "Enable MPI support" OFF)
opt(DTK "Enable DTK extensions" OFF)
opt(EXODUS "Enable EXODUS II extensions" OFF)

# set up version
set(NEMOSYS_MAJOR_VERSION 1)
set(NEMOSYS_MINOR_VERSION 2)
set(NEMOSYS_PATCH_VERSION 0)
set(NEMOSYS_EXTRA_VERSION "" CACHE STRING "NEMoSys extra version string")
set(NEMOSYS_VERSION "${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION}")
set(NEMOSYS_VERSION "${GMSH_VERSION}.${GMSH_PATCH_VERSION}${GMSH_EXTRA_VERSION}")
set(NEMOSYS_SHORT_LICENSE "GNU Lesser General Public License (LGPL)")

# finding libraries ###############################################################
IF (ENABLE_MPI)
  FIND_PACKAGE(MPI)
  IF(MPI_CXX_COMPILER)
    SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
  ENDIF()
ENDIF (ENABLE_MPI)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_LIBRARY(MADLIB_LIB MAdLib)
FIND_FILE(MADLIB_HDR MAdLib/MAdLib.h)
GET_FILENAME_COMPONENT(MADLIB_INCPATH ${MADLIB_HDR} PATH)

FIND_LIBRARY(GMSH_LIB Gmsh)
FIND_FILE(GMSH_HDR gmsh/Gmsh.h)
GET_FILENAME_COMPONENT(GMSH_INCPATH ${GMSH_HDR} PATH)
MESSAGE(STATUS ${GMSH_INCPATH})


FIND_LIBRARY(CGNS_LIB cgns)
FIND_FILE(CGNS_HDR cgnslib.h)
GET_FILENAME_COMPONENT(CGNS_INCPATH ${CGNS_HDR} PATH)

FIND_LIBRARY(ANN_LIB ann HINTS /usr/lib/)
FIND_FILE(ANN_HDR ANN.h HINTS /usr/include/ANN/)
GET_FILENAME_COMPONENT(ANN_INCPATH ${ANN_HDR} PATH)

FIND_LIBRARY(METIS_LIB metis)
FIND_FILE(METIS_HDR metis.h)
GET_FILENAME_COMPONENT(METIS_INCPATH ${METIS_HDR} PATH)

FIND_LIBRARY(NETGEN_LIB nglib HINTS ${CMAKE_SOURCE_DIR}/install/netgen/lib/)
FIND_FILE(NETGEN_HDR nglib.h HINTS ${CMAKE_SOURCE_DIR}/install/netgen/include/include)
GET_FILENAME_COMPONENT(NETGEN_INCPATH ${NETGEN_HDR} PATH)
MESSAGE(STATUS ${NETGEN_INCPATH})


IF (ENABLE_DTK)
  #FIND_LIBRARY(DTK_LIB classic_dtk)
  #FIND_FILE(DTK_HDR DTK_GeometryManager.hpp)
  #GET_FILENAME_COMPONENT(DTK_INCPATH ${DTK_HDR} PATH)
  #MESSAGE(${TRILINOS_PATH})
  FIND_PACKAGE(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})
  IF(Trilinos_FOUND)
     #MESSAGE("\nFound Trilinos!  Here are the details: ")
     #MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
     #MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
     #MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
     #MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
     #MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
     #MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
     #MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
     #MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
     #MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
     #MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
     #MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
     #MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
     #MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
     #MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
     #MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS =
     #  ${Trilinos_Fortran_COMPILER_FLAGS}")
     #MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
     #MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
     #MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
     #MESSAGE("End of Trilinos details\n")
  ELSE()
    MESSAGE(FATAL_ERROR "Could not find Trilinos!")
  ENDIF()
  #INCLUDE_DIRECTORIES(".")
  #INCLUDE_DIRECTORIES(${Trilinos_TPL_INCLUDE_DIRS})
  #INCLUDE_DIRECTORIES(${Trilinos_INCLUDE_DIRS})
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
INCLUDE_DIRECTORIES(AFTER ${ANN_INCPATH})
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

# Compiled objects output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set up RPATH
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Set up testing
OPTION(BUILD_TESTING "Build tests" OFF)
IF(BUILD_TESTING)
  INCLUDE(CTest)
  SET(BUILDNAME "${BUILDNAME}" CACHE STRING "Name of build on the dashboard")
  MARK_AS_ADVANCED(BUILDNAME)
ENDIF(BUILD_TESTING)


# Setting files   ###############################################################
# Install include files
FILE(GLOB H_FILES include/*.H)
FILE(GLOB h_FILES include/*.h)
FILE(GLOB hpp_FILES include/*.hpp)
SET(Nemosys_INC_FILES ${H_FILES} ${h_FILES} ${hpp_FILES})
INSTALL(FILES ${Nemosys_INC_FILES} DESTINATION include)

# Setting compile and link flags
SET(NEMOSYS_SRCS src/vtkAnalyzer.C src/cgnsAnalyzer.C src/meshPartitioner.C
                 src/cgnsWriter.C  src/gridTransfer.C src/meshPhys.C src/netgenInterface.C
                 src/Mesh/meshBase.C src/Mesh/vtkMesh.C src/MeshGeneration/meshGen.C 
                 src/Transfer/TransferBase.C  src/Transfer/FETransfer.C
                 src/SizeFieldGeneration/SizeFieldBase.C 
                 src/SizeFieldGeneration/GradSizeField.C 
                 src/SizeFieldGeneration/ValSizeField.C 
                 src/Refinement/Refine.C src/MeshQuality/MeshQuality.C 
                 src/Drivers/NemDriver.C src/Drivers/TransferDriver.C
                 src/Drivers/RefineDriver.C src/Drivers/MeshGenDriver.C 
                 src/Drivers/MeshQualityDriver.C src/StlToVtk.C src/Integration/Cubature.C)
SET(ROCSTAR_SRCS src/rocstarCgns.C)
SET(UTIL_SRCS utils/Nemosys.C utils/vol2planeTransfer.C 
    utils/cgns2msh.C utils/rocRemesh utils/rocSurfRemesh.C)
SET(TUTR_SRCS utils/writeStructured.C)
SET(R8_SRCS src/math/r8lib.cpp)
SET(SPHERE_SRCS src/spheres.cpp)
SET(INTERP_SRCS src/math/rbf_interp_nd.cpp src/math/rbfInterp.C src/math/baseInterp.C)
set(ALL_SRCS ${NEMOSYS_SRCS} ${UTIL_SRCS} ${TUTR_SRCS}
             ${R8_SRCS} ${INTERP_SRCS} ${ROCSTAR_SRCS} ${SPHERE_SRCS}
             testing/test_scripts/diff_vtu.cpp testing/test_scripts/diff_tol.cpp)
SET_SOURCE_FILES_PROPERTIES(${ALL_SRCS} PROPERTIES COMPILE_FLAGS "-fPIC -O3 ${MPI_CXX_COMPILE_FLAGS}" )
ADD_DEFINITIONS(-D_HAVE_GMSH_)

# Builiding libraries #########################################################

# build sphere lib
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(sphere SHARED ${SPHERE_SRCS})
ELSE()
  ADD_LIBRARY(sphere ${SPHERE_SRCS})
ENDIF()
INSTALL(TARGETS sphere LIBRARY DESTINATION lib)

# Build r8 and rbf library
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(r8 SHARED ${R8_SRCS})
ELSE()
  ADD_LIBRARY(r8 ${R8_SRCS})
ENDIF()
INSTALL(TARGETS r8 LIBRARY DESTINATION lib)

# build interp lib
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(interp SHARED ${INTERP_SRCS})
ELSE()
  ADD_LIBRARY(interp ${R8_SRCS})
ENDIF()
TARGET_LINK_LIBRARIES(interp r8)
INSTALL(TARGETS interp LIBRARY DESTINATION lib)

# build nemosys lib
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(Nemosys SHARED ${NEMOSYS_SRCS})
ELSE()
  ADD_LIBRARY(Nemosys ${NEMOSYS_SRCS})
  add_definitions( -DSTATIC_LINK )
ENDIF()
TARGET_LINK_LIBRARIES(Nemosys ${MADLIB_LIB} ${NETGEN_LIB} ${GMSH_LIB} ${VTK_LIBRARIES} ${CGNS_LIB} ${METIS_LIB} sphere interp)
INSTALL(TARGETS Nemosys LIBRARY DESTINATION lib)


# Rocstar support library
IF(NOT BUILD_STATIC)
  ADD_LIBRARY(Rocs SHARED ${ROCSTAR_SRCS})
ELSE()
  ADD_LIBRARY(Rocs ${ROCSTAR_SRCS})
ENDIF()
INSTALL(TARGETS Rocs LIBRARY DESTINATION lib)


# Builiding utilities #########################################################
# There are too many dependencies for gmsh this needs to be populated automatically
OPTION(BUILD_UTILS "BUILD_UTILS" OFF)

SET(CMAKE_CXX_LINK_EXECUTABLE_OLD ${CMAKE_CXX_LINK_EXECUTABLE})
SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lfltk_images -lfltk -lfltk_gl -lfltk_forms -lSM -lICE -lX11 -lXext -lm -ljpeg -lz -lpng -lGLU -lGL -llapack -lgmp -lblas -lcgns -lhdf5_cpp")

# Building Nemosys Core utility ###########################################
ADD_EXECUTABLE(nemosysRun utils/Nemosys.C)
if(VTK_LIBRARIES)
  TARGET_LINK_LIBRARIES(nemosysRun Nemosys)
else()
  TARGET_LINK_LIBRARIES(nemosysRun Nemosys vtkHybrid vtkWidgets)
endif()
INSTALL(TARGETS nemosysRun RUNTIME DESTINATION bin)

ADD_EXECUTABLE(nemTest utils/nemTest.C)
if(VTK_LIBRARIES)
  TARGET_LINK_LIBRARIES(nemTest Nemosys)
else()
  TARGET_LINK_LIBRARIES(nemTest Nemosys vtkHybrid vtkWidgets)
endif()
INSTALL(TARGETS nemTest RUNTIME DESTINATION bin)


IF(BUILD_UTILS)

  # xmlDump utility
  ADD_EXECUTABLE(xmlDump MACOSX_BUNDLE utils/xmlDump.C)
  if(VTK_LIBRARIES)
    TARGET_LINK_LIBRARIES(xmlDump Nemosys)
  else()
    TARGET_LINK_LIBRARIES(xmlDump Nemosys vtkHybrid vtkWidgets)
  ENDIF()
  INSTALL(TARGETS xmlDump RUNTIME DESTINATION bin)

  # cgns2gmsh utility
  ADD_EXECUTABLE(cgns2msh utils/cgns2msh.C)
  TARGET_LINK_LIBRARIES(cgns2msh Nemosys)
  INSTALL(TARGETS cgns2msh RUNTIME DESTINATION bin)

  # Rocstart volume remesh utility
  ADD_EXECUTABLE(rocRemesh utils/rocRemesh.C)
  TARGET_LINK_LIBRARIES(rocRemesh Nemosys)
  INSTALL(TARGETS rocRemesh RUNTIME DESTINATION bin)

  # Rocstar surface remesh utility
  ADD_EXECUTABLE(rocSurfRemesh utils/rocSurfRemesh.C)
  TARGET_LINK_LIBRARIES(rocSurfRemesh Nemosys Rocs)
  INSTALL(TARGETS rocSurfRemesh RUNTIME DESTINATION bin)

  # Rocstar surface remesh utility
  ADD_EXECUTABLE(grid2gridTransfer utils/grid2gridTransfer.C)
  TARGET_LINK_LIBRARIES(grid2gridTransfer Nemosys Rocs)
  INSTALL(TARGETS grid2gridTransfer RUNTIME DESTINATION bin)
  
  #vol2planeTransfer utility
  ADD_EXECUTABLE(vol2planeTransfer utils/vol2planeTransfer.C)
  if(VTK_LIBRARIES)
    TARGET_LINK_LIBRARIES(vol2planeTransfer Nemosys)
  else()
    TARGET_LINK_LIBRARIES(vol2planeTransfer Nemosys vtkHybrid vtkWidgets)
  ENDIF()
  INSTALL(TARGETS vol2planeTransfer RUNTIME DESTINATION bin)

 # Building tutorials ###########################################################
 SET(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE_OLD} -lhdf5_cpp")

  # meshTransfer tutorial
  ADD_EXECUTABLE(mshTransfer tutorials/mshTransfer.C)
  TARGET_LINK_LIBRARIES(mshTransfer Nemosys)
  INSTALL(TARGETS mshTransfer RUNTIME DESTINATION bin)

  # DTK tutorial
  IF (ENABLE_DTK)
    ADD_EXECUTABLE(dtkTutorial tutorials/dtkTutorial.C)
    SET_TARGET_PROPERTIES(dtkTutorial PROPERTIES COMPILE_FLAGS ${dtk_headers} LINK_FLAGS ${dtk_libs})
    TARGET_LINK_LIBRARIES(dtkTutorial Nemosys Rocs
                                             ${DTK_LIB} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
    INSTALL(TARGETS dtkTutorial RUNTIME DESTINATION bin)
  ENDIF()

ENDIF(BUILD_UTILS) 



# Building converters ##########################################################
IF (ENABLE_EXODUS)
  ADD_SUBDIRECTORY(converter)
ENDIF ()

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
  INSTALL( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs DESTINATION share/Nemosys )
ENDIF()

# Perform testing ############################################################
# adding gtest subdir with its own cmake file
add_subdirectory(${CMAKE_SOURCE_DIR}/testing/lib/gtest-1.6.0)
# include gtest library. gtest_SOURCE_DIR is availble from subdir addition
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}/include/gtest ${gtest_SOURCE_DIR})

SET(CMAKE_CXX_LINK_EXECUTABLE ${CMAKE_CXX_LINK_EXECUTABLE_OLD})
SET(DIFFTOL ${CMAKE_SOURCE_DIR}/testing/test_scripts/diff_tol)
SET(DIFFTOLSOURCE ${CMAKE_SOURCE_DIR}/testing/test_scripts/diff_tol.cpp)
SET(DIFFVTU ${CMAKE_SOURCE_DIR}/testing/test_scripts/diff_vtu)
SET(DIFFVTUSOURCE ${CMAKE_SOURCE_DIR}/testing/test_scripts/diff_vtu.cpp)
SET(OLD_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/testing/test_scripts)
ADD_EXECUTABLE(diff_tol ${DIFFTOLSOURCE})
ADD_EXECUTABLE(diff_vtu ${DIFFVTUSOURCE})
ADD_EXECUTABLE(runCubatureInterpTest ${CMAKE_SOURCE_DIR}/testing/test_scripts/testCubature.C)
ADD_EXECUTABLE(runConversionTest ${CMAKE_SOURCE_DIR}/testing/test_scripts/testConversion.C)
TARGET_LINK_LIBRARIES(runCubatureInterpTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runConversionTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(diff_vtu Nemosys)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OLD_RUNTIME_OUTPUT_DIRECTORY})
ADD_SUBDIRECTORY(testing)
INSTALL(FILES ${DIFFTOL} ${DIFFVTU} PERMISSIONS OWNER_EXECUTE OWNER_READ OWNER_WRITE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE DESTINATION ${CMAKE_SOURCE_DIR}/testing/test_scripts)
INSTALL(TARGETS runCubatureInterpTest RUNTIME DESTINATION ${CMAKE_SOURCE_DIR}/testing/test_scripts)
