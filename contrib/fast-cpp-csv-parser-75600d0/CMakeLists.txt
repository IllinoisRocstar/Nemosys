# https://github.com/ben-strasser/fast-cpp-csv-parser
cmake_minimum_required(VERSION 3.9)
project(fast-cpp-csv-parser)

add_library(fast-cpp-csv-parser INTERFACE)

include(GNUInstallDirs)

target_include_directories(fast-cpp-csv-parser INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/fast-cpp-csv-parser>)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
if(Threads_FOUND)
  target_link_libraries(fast-cpp-csv-parser INTERFACE Threads::Threads)
else()
  target_compile_definitions(fast-cpp-csv-parser INTERFACE CSV_IO_NO_THREAD)
endif()

install(FILES csv.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fast-cpp-csv-parser)

install(TARGETS fast-cpp-csv-parser EXPORT Nemosys)
