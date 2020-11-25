include(CMakeFindDependencyMacro)
find_dependency(jsoncons)
find_dependency(Omega_h)

include("${CMAKE_CURRENT_LIST_DIR}/Nemosys.cmake")
