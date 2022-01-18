# https://github.com/akuukka/quickhull
cmake_minimum_required(VERSION 3.9)
project(quickhull)

add_library(quickhull INTERFACE)

include(GNUInstallDirs)

target_include_directories(quickhull INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/quickhull>)

install(DIRECTORY .
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/quickhull
    FILES_MATCHING PATTERN "*.hpp"
                   PATTERN "*.cpp"
                   PATTERN "Tests" EXCLUDE)

install(TARGETS quickhull EXPORT Nemosys)
