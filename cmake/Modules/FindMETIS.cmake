find_library(METIS_LIBRARY metis)
find_path(METIS_INCLUDE_DIR NAMES metis.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
    REQUIRED_VARS METIS_INCLUDE_DIR METIS_LIBRARY)

set(HAVE_METIS TRUE)
add_library(METIS::METIS UNKNOWN IMPORTED)
set_target_properties(METIS::METIS PROPERTIES
    IMPORTED_LOCATION "${METIS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}")
