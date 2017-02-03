cmake_minimum_required (VERSION 2.8.11)
project (RBFINTERP)

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single directory for all executables.")
set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib CACHE PATH "Single directory for all libraries and archives.")
mark_as_advanced (LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# include files
include_directories(include)

# build r8 library
add_library(r8 SHARED src/r8lib.cpp)

# build rbf library
add_library(rbf SHARED src/rbf_interp_nd.cpp src/rbfInterp.C)

# build test executables
add_executable(test test/rbf_interp_nd_prb.cpp)
target_link_libraries(test r8 rbf)
add_executable(rbf_class_test test/rbfClassTest.C)
target_link_libraries(rbf_class_test r8 rbf)

# instllation directives
install (TARGETS r8 rbf DESTINATION lib)
install (TARGETS test DESTINATION bin)
