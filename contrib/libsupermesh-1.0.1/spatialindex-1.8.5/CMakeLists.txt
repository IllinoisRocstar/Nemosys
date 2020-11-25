cmake_minimum_required(VERSION 3.1)
project(spatialindex)

# internal cmake settings #############################################
set(CMAKE_COLOR_MAKEFILE ON)
mark_as_advanced(CMAKE_VERBOSE_MAKEFILE)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

if (CMAKE_MAJOR_VERSION GREATER 2)
    cmake_policy(SET CMP0042 OLD) # osx rpath
endif()

# project general settings ############################################
SET(SIDX_VERSION_MAJOR "1")
SET(SIDX_VERSION_MINOR "8")
SET(SIDX_VERSION_PATCH "5")
set(SIDX_VERSION_STRING "${SIDX_VERSION_MAJOR}.${SIDX_VERSION_MINOR}.${SIDX_VERSION_PATCH}")

# TODO: this may need to be removed depending on how windows compiles
if(WIN32)
  if (MSVC)
    if (MSVC80 OR MSVC90 OR MSVC10 OR MSVC11 OR MSVC12 OR MSVC13)
        if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
            set( SIDX_LIB_NAME "spatialindex-64" )
            set( SIDX_C_LIB_NAME "spatialindex_c-64" )
        else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
            set( SIDX_LIB_NAME "spatialindex-32"  )
            set( SIDX_C_LIB_NAME "spatialindex_c-32"  )
        endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    endif()
  endif()
endif()


include(CheckFunctionExists)
check_function_exists(srand48 HAVE_SRAND48)
check_function_exists(gettimeofday HAVE_GETTIMEOFDAY)
check_function_exists(memset HAVE_MEMSET)
check_function_exists(memcpy HAVE_MEMCPY)
check_function_exists(bcopy HAVE_BCOPY)

include(CheckIncludeFiles)
check_include_files(pthread.h HAVE_PTHREAD_H)

set(SIDX_COMPILE_DEFS "")

if(HAVE_SRAND48)
  list(APPEND SIDX_COMPILE_DEFS "-DHAVE_SRAND48")
endif()

if (HAVE_GETTIMEOFDAY)
  list(APPEND SIDX_COMPILE_DEFS "-DHAVE_GETTIMEOFDAY")
endif()

if (HAVE_BZERO)
  list(APPEND SIDX_COMPILE_DEFS "-DHAVE_BZERO")
endif()

if (HAVE_MEMSET)
  list(APPEND SIDX_COMPILE_DEFS "-DHAVE_MEMSET")
endif()

if (HAVE_MEMCPY)
  list(APPEND SIDX_COMPILE_DEFS "-DHAVE_MEMCPY")
endif()

if (HAVE_BCOPY)
  list(APPEND SIDX_COMPILE_DEFS "-DHAVE_BCOPY")
endif()

if (HAVE_PTHREAD_H)
  list(APPEND SIDX_COMPILE_DEFS "-DHAVE_PTHREAD_H")
endif()


# OLD

#set(SIDX_COMPILE_DEFS "")
#if (HAVE_SRAND48)
#  list(APPEND SIDX_COMPILE_DEFS HAVE_SRAND48)
#endif()
#
#if (HAVE_GETTIMEOFDAY)
#  list(APPEND SIDX_COMPILE_DEFS HAVE_GETTIMEOFDAY)
#endif()
#
#if (HAVE_BZERO)
#  list(APPEND SIDX_COMPILE_DEFS HAVE_BZERO)
#endif()
#
#if (HAVE_MEMSET)
#  list(APPEND SIDX_COMPILE_DEFS HAVE_MEMSET)
#endif()
#
#if (HAVE_MEMCPY)
#  list(APPEND SIDX_COMPILE_DEFS HAVE_MEMCPY)
#endif()
#
#if (HAVE_BCOPY)
#  list(APPEND SIDX_COMPILE_DEFS HAVE_BCOPY)
#endif()
#
#if (HAVE_PTHREAD_H)
#  list(APPEND SIDX_COMPILE_DEFS HAVE_PTHREAD_H)
#endif()
#

# TODO: this may need to be removed depending on how windows compiles
if(WIN32)
  if (MSVC)
    if (MSVC80 OR MSVC90 OR MSVC10 OR MSVC11 OR MSVC12 OR MSVC13)
      add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
      add_definitions(-D_CRT_SECURE_NO_WARNINGS)
      add_definitions(-D_CRT_NONSTDC_NO_WARNING)
      add_definitions(-D_SCL_SECURE_NO_WARNINGS)
      add_definitions(-DSPATIALINDEX_CREATE_DLL=1)
      add_definitions(-DNOMINMAX)

      set(SIDX_COMPILER_MSVC 1)
    endif()

    if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
      string(REGEX REPLACE "/W[0-4]" "/W4"
        CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    endif()
  endif()

  set(CMAKE_INCLUDE_PATH c:/osgeo4w/include;$ENV{CMAKE_INCLUDE_PATH})
  set(CMAKE_LIBRARY_PATH c:/osgeo4w/lib;$ENV{CMAKE_LIBRARY_PATH})

  set(SIDX_OS_WIN32 1)

else() # not WIN32

  # Recommended C++ compilation flags -Weffc++
  set(SIDX_COMMON_CXX_FLAGS
    "-pedantic -Wall -Wpointer-arith -Wcast-align -Wcast-qual  -Wredundant-decls -Wno-long-long")

  if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIDX_COMMON_CXX_FLAGS}")
    if (CMAKE_COMPILER_IS_GNUCXX)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98 -ansi")
    endif()
    set(SIDX_COMPILER_GCC 1)
  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER}" MATCHES "clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SIDX_COMMON_CXX_FLAGS}")
    set(SIDX_COMPILER_CLANG 1)
  endif()
endif(WIN32)

if (APPLE)
  set(SO_EXT dylib)
  set(CMAKE_FIND_FRAMEWORK "LAST")
elseif(WIN32)
  set(SO_EXT dll)
else()
  set(SO_EXT so)
endif(APPLE)

# installation path settings
if(WIN32)
  set(DEFAULT_LIB_SUBDIR lib)
  set(DEFAULT_DATA_SUBDIR .)
  set(DEFAULT_INCLUDE_SUBDIR include)

  if (MSVC)
    set(DEFAULT_BIN_SUBDIR bin)
  else()
    set(DEFAULT_BIN_SUBDIR .)
  endif()
else()
  # Common locations for Unix and Mac OS X
  set(DEFAULT_BIN_SUBDIR bin)
  set(DEFAULT_LIB_SUBDIR lib${LIB_SUFFIX})
  set(DEFAULT_DATA_SUBDIR share/spatialindex)
  set(DEFAULT_INCLUDE_SUBDIR include)
endif()


# building targets ################################################################
add_subdirectory(src)
# tests turned off
#add_subdirectory(test)
