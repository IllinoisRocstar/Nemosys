set(Geom_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

add_library(Geom ${Geom_SRC_DIR}/GEOMUtils/GEOMUtils.cxx)
target_link_libraries(Geom PUBLIC Kernel TKBRep TKShHealing TKMesh TKGeomAlgo TKTopAlgo TKG2d TKG3d TKV3d TKGeomBase TKBO)
target_include_directories(Geom PUBLIC
  $<BUILD_INTERFACE:${Geom_SRC_DIR}/GEOMUtils>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/Geom>
  ${OpenCASCADE_INCLUDE_DIR}
)

if(WIN32)
    target_compile_definitions(Geom PRIVATE UTILS_EXPORTS)
    target_compile_definitions(Geom PRIVATE UNICODE _UNICODE)
endif()

install(FILES ${Geom_SRC_DIR}/GEOMUtils/GEOMUtils.hxx DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/Geom)

set(Geom_LIBRARIES Geom)

install(TARGETS ${Geom_LIBRARIES} EXPORT SMESH-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
