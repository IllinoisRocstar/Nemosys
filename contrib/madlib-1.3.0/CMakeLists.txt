cmake_minimum_required(VERSION 2.8)
project(MAdLib)

message(STATUS "MAdLib Base")

add_subdirectory(Adapt)
add_subdirectory(Benchmarks/checkMesh)
add_subdirectory(Benchmarks/meshInfo)
add_subdirectory(Benchmarks/moveIt)
add_subdirectory(Benchmarks/optimize)
add_subdirectory(Common)
add_subdirectory(Contrib/ANN)
add_subdirectory(Geo)
add_subdirectory(Mesh)