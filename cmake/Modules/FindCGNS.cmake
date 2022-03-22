find_library(CGNS_LIBRARY cgns)
find_path(CGNS_INCLUDE_DIR NAMES cgnslib.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CGNS
    REQUIRED_VARS CGNS_INCLUDE_DIR CGNS_LIBRARY)

set(HAVE_CGNS TRUE)
add_library(CGNS::CGNS UNKNOWN IMPORTED)
set_target_properties(CGNS::CGNS PROPERTIES
    IMPORTED_LOCATION "${CGNS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CGNS_INCLUDE_DIR}")
