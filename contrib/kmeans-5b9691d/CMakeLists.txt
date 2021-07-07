# https://github.com/felixduvallet/kmeans
cmake_minimum_required(VERSION 3.9)
project(kmeans)

set(KMEANS_PUBLIC_HEADERS
    src/point.h
    src/kmeans.h)

set(KMEANS_SOURCE_FILES
    src/point.cpp
    src/kmeans.cpp
    ${KMEANS_PUBLIC_HEADERS})

add_library(kmeans_lib ${KMEANS_SOURCE_FILES})

set_target_properties(kmeans_lib PROPERTIES CXX_STANDARD 11)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  set_target_properties(kmeans_lib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

include(GNUInstallDirs)

target_include_directories(kmeans_lib
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/kmeans>)

install(FILES ${KMEANS_PUBLIC_HEADERS}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/kmeans)

install(TARGETS kmeans_lib
    EXPORT Nemosys
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
