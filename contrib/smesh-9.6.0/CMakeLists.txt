cmake_minimum_required(VERSION 3.9)
project(SMESH C CXX)

# --------------------------------------------------------------------------- #
# OPTIONS
# --------------------------------------------------------------------------- #
option(ENABLE_NETGEN "Enable Netgen" ON)
option(ENABLE_GMSH "Enable Gmsh" ON)
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE PATH "Installation directory")
set(CMAKE_INSTALL_LIBDIR lib CACHE PATH "Output directory for libraries")


# --------------------------------------------------------------------------- #
# SETTINGS
# --------------------------------------------------------------------------- #
set(SMESH_VERSION_MAJOR 9)
set(SMESH_VERSION_MINOR 6)
set(SMESH_VERSION_PATCH 0)
set(SMESH_VERSION_TWEAK 3)

# Build shared libraries
set(BUILD_SHARED_LIBS TRUE)

# Force C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Warning and errors
if(MSVC)
    if(POLICY CMP0092)
        # CMP0092 introduced in CMake 3.15. Removes /W3 as default warning flag
        # for MSVC-based compilers.
        cmake_policy(SET CMP0092 NEW)
        add_compile_options(/W1)
    else()
        string(REGEX REPLACE "/W3" "/W1" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    endif()
    add_compile_options(/wd4715 /wd4091 /wd4503 /wd4805 /wd4267 /wd4477 /wd4138 /wd4273 /wd4311 /wd4806)
else()
    add_compile_options(-Wno-deprecated)
endif()

# Definitions
if(UNIX)
  # Same settings are used for both MacOSX and Unix/Linux
  add_definitions(-DHAVE_LIMITS_H -DCSFDB -DLIN -DOCC_CONVERT_SIGNALS)
else(UNIX)
  if(WIN32)
      if(MSVC)
          add_definitions(-DWNT -DWIN32 -D_WINDOWS -DCSFDB -DUSE_CLOCK -DMSDOS -DNO_ONEXIT -DNO_My_ctype -DNO_ISATTY -DNO_FPINIT)
      else(MSVC)
          add_definitions(-DWNT -DWIN32 -D_WINDOWS -DCSFDB)
      endif(MSVC)
  else(WIN32)
      message(FATAL_ERROR "Unknown platform")
  endif(WIN32)
endif(UNIX)

if(CMAKE_SIZEOF_VOID_P STREQUAL 8)
    add_definitions(-D_OCC64)
endif(CMAKE_SIZEOF_VOID_P STREQUAL 8)


# --------------------------------------------------------------------------- #
# OpenCASCADE
# --------------------------------------------------------------------------- #
message(STATUS "Searching for OpenCASCADE...")
find_package(OpenCASCADE REQUIRED)


# --------------------------------------------------------------------------- #
# VTK
# --------------------------------------------------------------------------- #
message(STATUS "Searching for VTK...")
find_package(VTK REQUIRED)
if(VTK_VERSION VERSION_GREATER_EQUAL 8.90)
  find_package(VTK REQUIRED COMPONENTS CommonCore CommonDataModel FiltersVerdict)
else()
  # Would add as alias as well, but need to support older CMake
  find_package(VTK REQUIRED COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersVerdict)
endif()


# --------------------------------------------------------------------------- #
# BOOST
# --------------------------------------------------------------------------- #
message(STATUS "Searching for Boost...")
if(WIN32)
  set(Boost_USE_STATIC_LIBS ON)
endif()
set(BOOST_USE_MULTITHREADED ON)
find_package(Boost REQUIRED COMPONENTS filesystem thread serialization)

if (MSVC)
  # set postfix for debug libs
  if(NOT CMAKE_DEBUG_POSTFIX)
      set(CMAKE_DEBUG_POSTFIX d)
  endif()
endif()


# --------------------------------------------------------------------------- #
# SMESH
# --------------------------------------------------------------------------- #

include(GNUInstallDirs)

# Kernel
add_subdirectory(Kernel)

# Geom
add_subdirectory(Geom)

# SMESH
add_subdirectory(SMESH)

# NETGENPlugin
if(ENABLE_NETGEN)
  find_package(Netgen REQUIRED)
  add_subdirectory(NETGENPlugin)
endif(ENABLE_NETGEN)

# GMSHPlugin
if(ENABLE_GMSH)
  find_package(Gmsh REQUIRED)
  add_subdirectory(GMSHPlugin)
endif(ENABLE_GMSH)

# --------------------------------------------------------------------------- #
# Install
# --------------------------------------------------------------------------- #
install(EXPORT SMESH-targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smesh NAMESPACE SMESH::)


# --------------------------------------------------------------------------- #
# Configuration files
# --------------------------------------------------------------------------- #
include(CMakePackageConfigHelpers)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/SMESHConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/SMESHConfig.cmake @ONLY)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/SMESHConfigVersion.cmake
  VERSION ${SMESH_VERSION_MAJOR}.${SMESH_VERSION_MINOR}.${SMESH_VERSION_PATCH}
  COMPATIBILITY ExactVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SMESHConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/SMESHConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/smesh)
