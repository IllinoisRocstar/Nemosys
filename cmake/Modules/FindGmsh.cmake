# This module defines
#  Gmsh_INCLUDE_DIRS: path to gmsh.h (public header)
#  Gmsh_PRIVATE_INCLUDE_DIRS: path to GmshVersion.h (private header); generally ${Gmsh_INCLUDE_DIR}/gmsh
#  Gmsh_LIBRARIES: the Gmsh library
#  Gmsh_FOUND: if false, do not try to use Gmsh
#  Gmsh_VERSION: The found version of Gmsh
#  Gmsh::Gmsh target

find_path(Gmsh_INCLUDE_DIR NAMES gmsh.h PATH_SUFFIXES gmsh/include)
mark_as_advanced(Gmsh_INCLUDE_DIR)

find_path(Gmsh_PRIVATE_INCLUDE_DIR GmshVersion.h HINTS ${Gmsh_INCLUDE_DIR}/gmsh)
mark_as_advanced(Gmsh_PRIVATE_INCLUDE_DIR)

find_library(Gmsh_LIBRARY NAMES gmsh PATH_SUFFIXES gmsh/lib)
mark_as_advanced(Gmsh_LIBRARY)

set(Gmsh_VERSION Gmsh_VERSION-NOTFOUND)
if(Gmsh_PRIVATE_INCLUDE_DIR)
  if(EXISTS "${Gmsh_PRIVATE_INCLUDE_DIR}/GmshVersion.h")
    file(STRINGS "${Gmsh_PRIVATE_INCLUDE_DIR}/GmshVersion.h" _gmsh_version REGEX "GMSH_.*_VERSION")
    string(REGEX REPLACE ".*GMSH_MAJOR_VERSION *\([0-9]*\).*" "\\1" _gmsh_major "${_gmsh_version}")
    string(REGEX REPLACE ".*GMSH_MINOR_VERSION *\([0-9]*\).*" "\\1" _gmsh_minor "${_gmsh_version}")
    string(REGEX REPLACE ".*GMSH_PATCH_VERSION *\([0-9]*\).*" "\\1" _gmsh_patch "${_gmsh_version}")
    unset(_gmsh_version)
    if(NOT _gmsh_major STREQUAL "" AND NOT _gmsh_minor STREQUAL "" AND NOT _gmsh_patch STREQUAL "")
      set(Gmsh_VERSION "${_gmsh_major}.${_gmsh_minor}.${_gmsh_patch}")
    endif()
    unset(_gmsh_major)
    unset(_gmsh_minor)
    unset(_gmsh_patch)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Gmsh
    REQUIRED_VARS Gmsh_LIBRARY Gmsh_INCLUDE_DIR Gmsh_PRIVATE_INCLUDE_DIR
    VERSION_VAR Gmsh_VERSION)

if(Gmsh_FOUND)
  set(Gmsh_LIBRARIES "${Gmsh_LIBRARY}")
  set(Gmsh_INCLUDE_DIRS "${Gmsh_INCLUDE_DIR}")
  set(Gmsh_PRIVATE_INCLUDE_DIRS "${Gmsh_PRIVATE_INCLUDE_DIR}")
  if(NOT TARGET Gmsh::Gmsh)
    add_library(Gmsh::Gmsh UNKNOWN IMPORTED)
    set_target_properties(Gmsh::Gmsh PROPERTIES
        IMPORTED_LOCATION "${Gmsh_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${Gmsh_INCLUDE_DIR}"
        INTERFACE_COMPILE_DEFINITIONS $<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>)
  endif()
endif()
