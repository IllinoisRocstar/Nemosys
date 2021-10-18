set(ANN_MAJOR_VERSION 1)
set(ANN_MINOR_VERSION 1)
set(ANN_PATCH_VERSION 2)
#set(ANN_EXTRA_VERSION "")
set(ANN_VERSION "${ANN_MAJOR_VERSION}.${ANN_MINOR_VERSION}")
set(ANN_VERSION "${ANN_VERSION}.${ANN_PATCH_VERSION}${ANN_EXTRA_VERSION}")

# Set ANN include files
set(ANN_INC_FILES
    include/ANN/ANN.h
#    include/ANN/ANNperf.h
#    include/ANN/ANNx.h
)

# Set ANN source files
set(ANN_SRCS
    src/ANN.cpp
    src/bd_fix_rad_search.cpp
    src/bd_pr_search.cpp
    src/bd_search.cpp
    src/bd_tree.cpp
    src/brute.cpp
    src/kd_dump.cpp
    src/kd_fix_rad_search.cpp
    src/kd_pr_search.cpp
    src/kd_search.cpp
    src/kd_split.cpp
    src/kd_tree.cpp
    src/kd_util.cpp
    src/perf.cpp
)

# Build ANN library
add_library(ANN ${ANN_SRCS})

# Set ANN include directory
target_include_directories(ANN
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

set_target_properties(ANN PROPERTIES DEFINE_SYMBOL DLL_EXPORTS)

# Set soname version
set_target_properties(ANN PROPERTIES VERSION ${ANN_VERSION}
                                     SOVERSION ${ANN_MAJOR_VERSION})

# Install ANN headers
install(FILES ${ANN_INC_FILES}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ANN)

# Install ANN library
install(TARGETS ANN
    EXPORT Nemosys
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
