message(STATUS "Benchmarks/meshInfo")

add_executable(meshInfo main.cc)
target_link_libraries(meshInfo PRIVATE MAdLib)
