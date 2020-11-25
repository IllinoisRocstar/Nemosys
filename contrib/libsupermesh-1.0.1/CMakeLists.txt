cmake_minimum_required(VERSION 3.1)

project(LIBSUPERMESH C CXX Fortran)
#enable_language(Fortran)

# Set project metadata #########################################################
set(LIBSUPERMESH_VERSION_MAJOR 1)
set(LIBSUPERMESH_VERSION_MINOR 0)
set(LIBSUPERMESH_VERSION_PATCH 1)

# Set project options #########################################################
option(LIBSUPERMESH_ENABLE_JUDY "Enable use of the Judy library" OFF)
option(LIBSUPERMESH_ENABLE_TIMERS "Enable internal timers" OFF)
option(LIBSUPERMESH_OVERLAP_COMPUTE_COMMS "Overlap computation and communication. Efficiency of this depends upon the specific MPI implementation." OFF)
option(LIBSUPERMESH_DOUBLE_PRECISION "Build with double precision Fortran reals" ON)
option(LIBSUPERMESH_AUTO_COMPILER_FLAGS "Choose compiler flags automatically. This overrides any other compiler flags defined via CMake variables." ON)

# Set project dependencies #########################################################
# TODO: MPI is assumed to be optional
#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})
#if(MPI_LIBRARY)
#  set(MPI_LIBRARIES ${MPI_LIBRARY})
#endif(MPI_LIBRARY)
#if(MPI_EXTRA_LIBRARY)
#  set(MPI_LIBRARIES ${MPI_LIBRARIES} "${MPI_EXTRA_LIBRARY}")
#endif(MPI_EXTRA_LIBRARY)
#

if(LIBSUPERMESH_ENABLE_JUDY)
  find_library(JUDY_LIBRARY NAMES Judy)
  if(JUDY_LIBRARY STREQUAL JUDY_LIBRARY-NOTFOUND)
    message(FATAL_ERROR "Judy library not found")
  endif()
  list(APPEND link_libraries ${JUDY_LIBRARY})
  find_path(JUDY_INCLUDE_DIR NAMES Judy.h)
  if(JUDY_INCLUDE_DIR STREQUAL JUDY_INCLUDE_DIR-NOTFOUND)
    message(FATAL_ERROR "Judy header file not found")
  endif()
  include_directories(${JUDY_INCLUDE_DIR})
endif()

find_package(Backtrace)
if(Backtrace_FOUND)
  set(link_libraries ${link_libraries} ${Backtrace_LIBRARIES})
  include_directories(${Backtrace_INCLUDE_DIRS})
endif()

# Set project sources ##############################################################
file(GLOB LIBSMSRC
    ${LIBSUPERMESH_SOURCE_DIR}/src/*.F90
    ${LIBSUPERMESH_SOURCE_DIR}/src/*.c
    ${LIBSUPERMESH_SOURCE_DIR}/src/*.cpp)

# Set project compile flags #########################################################
set(CMAKE_Fortran_MODULE_DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/include)


if(LIBSUPERMESH_AUTO_COMPILER_FLAGS)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-ffree-line-length-none -pipe -std=f2008")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -Wall -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-integer=-66666 -finit-real=nan -fimplicit-none")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(CMAKE_Fortran_FLAGS "-s -fpic")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  endif()

  if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "-pipe -std=c90")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    endif()
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g -Wall")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG")
  elseif(CMAKE_C_COMPILER_ID MATCHES "Cray")
    set(CMAKE_C_FLAGS "-fpic")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_C_FLAGS_RELEASE "-O2")
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-pipe -std=c++03")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    endif()
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ftree-vectorize -DNDEBUG")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    set(CMAKE_CXX_FLAGS "-fpic")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
  endif()
else()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -std=f2008")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -s -fpic")
  endif()

  if(CMAKE_C_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c90")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    endif()
  elseif(CMAKE_C_COMPILER_ID MATCHES "Cray")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic")
  endif()

  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++03")
    if(BUILD_SHARED_LIBS)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    endif()
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic")
  endif()
endif()

message(STATUS "CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE=${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG=${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")

configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/config/libsupermesh.pc.in"
  "${LIBSUPERMESH_BINARY_DIR}/config/libsupermesh.pc"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh.h.in"
  "${LIBSUPERMESH_BINARY_DIR}/include/libsupermesh.h"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh-c.h"
  "${LIBSUPERMESH_BINARY_DIR}/include/libsupermesh-c.h"
  )
configure_file (
  "${LIBSUPERMESH_SOURCE_DIR}/include/libsupermesh_configuration.h.in"
  "${LIBSUPERMESH_BINARY_DIR}/include_local/libsupermesh_configuration.h"
  )

# Build dependencies  ########################################################
add_subdirectory(spatialindex-1.8.5)

# Build main library #########################################################
add_library(supermesh ${LIBSMSRC})

target_link_libraries(supermesh
    PUBLIC
        ${MPI_LIBRARIES}
        spatialindex)

target_include_directories(supermesh
    PUBLIC
        $<INSTALL_INTERFACE:include/supermesh>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include_local>
        ${MPI_Fortran_INCLUDE_PATH}
        ${MPI_INCLUDE_PATH} )

set_target_properties(supermesh
    PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}" )

if(CMAKE_BUILD_TYPE MATCHES "debug")
  target_compile_definitions(supermesh
      PUBLIC
      LIBSUPERMESH_DEBUG )
endif()

# turnning all tests off
#set(test_main ${LIBSUPERMESH_SOURCE_DIR}/src/tests/test_main.cpp)
#file(GLOB unittest_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/*.F90)
#set(unittests "")
#foreach(unittest_file ${unittest_files})
#  get_filename_component(unittest ${unittest_file} NAME_WE)
#  add_executable(${unittest} ${unittest_file} ${test_main})
#  list(APPEND unittests ${unittest})
#  set_property(TARGET ${unittest} PROPERTY COMPILE_DEFINITIONS "TESTNAME=${unittest}")
#  if(unittest MATCHES parallel)
#    add_test(${unittest} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${unittest} ${MPIEXEC_POSTFLAGS})
#  else()
#    add_test(${unittest} ${unittest})
#  endif()
#  set_tests_properties(${unittest} PROPERTIES FAIL_REGULAR_EXPRESSION "Fail:")
#endforeach()
#set(test_link_libraries supermesh ${link_libraries} ${MPI_CXX_LIBRARIES})
#target_link_libraries(supermesh ${link_libraries})
#foreach(test ${unittests})
#  target_link_libraries(${test} ${test_link_libraries})
#endforeach()

#set_target_properties(supermesh PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#foreach(unittest ${unittests})
#  set_target_properties(${unittest} PROPERTIES LINK_FLAGS "${MPI_Fortran_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS}")
#endforeach()
#
#option(ENABLE_DOCS "Enable building of documentation" OFF)
#if(ENABLE_DOCS)
#  find_package(ImageMagick COMPONENTS convert REQUIRED)
#  find_package(LATEX)
#  if(NOT LATEX_FOUND OR NOT LATEX_PDFLATEX_FOUND OR NOT LATEX_BIBTEX_FOUND)
#    message(SEND_ERROR "Required LaTeX components not found")
#  endif()
#
#  include(${LIBSUPERMESH_SOURCE_DIR}/cmake/UseLATEX.cmake)
#  set(LATEX_OUTPUT_PATH ${LIBSUPERMESH_BINARY_DIR}/doc)
#  add_latex_document(doc/manual.tex INPUTS doc/version.tex
#    CONFIGURE doc/version.tex IMAGES doc/hex_ordering.pdf TARGET_NAME doc
#    FORCE_PDF EXCLUDE_FROM_ALL BIBFILES doc/bibliography.bib)
#endif()

#install(TARGETS supermesh DESTINATION lib)
#install(DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
#install(FILES ${LIBSUPERMESH_BINARY_DIR}/config/libsupermesh.pc DESTINATION lib/pkgconfig)

#file(GLOB data_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.ele
#                     ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.node)
#file(COPY ${data_files} DESTINATION ${LIBSUPERMESH_BINARY_DIR}/data)
#file(GLOB data_files ${LIBSUPERMESH_SOURCE_DIR}/src/tests/data/*.tar.bz2)
#foreach(data_file ${data_files})
#  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xjf ${data_file} WORKING_DIRECTORY ${LIBSUPERMESH_BINARY_DIR}/data)
#endforeach()

# Install MAdLib library
install(TARGETS supermesh
    EXPORT Nemosys
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
