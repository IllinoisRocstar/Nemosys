message(STATUS "Benchmarks/checkMesh")

add_executable(checkMesh main.cc)
target_link_libraries(checkMesh PRIVATE MAdLib)
