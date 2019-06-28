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

# Generate export header #######################################################
include(GenerateExportHeader)
generate_export_header(ANN)
list(APPEND ${ANN_INC_FILES} ann_export.h)
file(COPY ${CMAKE_CURRENT_BINARY_DIR}/ann_export.h
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/include/ANN)

# Install ANN headers
install(FILES ${ANN_INC_FILES}
    DESTINATION include/ANN)

# Install ANN library
install(TARGETS ANN
    EXPORT Nemosys
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
