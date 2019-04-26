message(STATUS "Benchmarks/moveIt")

add_executable(moveIt main.cc)
target_link_libraries(moveIt PRIVATE MAdLib)
