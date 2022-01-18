# GMSHPlugin
set(GMSHPlugin_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/GMSHPlugin)

# Get sources and remove unsupported ones
file(GLOB GMSHPlugin_SRC ${GMSHPlugin_SRC_DIR}/*.cxx)
add_library(GMSHPlugin ${GMSHPlugin_SRC})
target_link_libraries(GMSHPlugin PUBLIC Kernel SMESH StdMeshers SMDS Gmsh::Gmsh)
target_link_libraries(GMSHPlugin PRIVATE Controls SMESHDS TKernel TKBRep TKMath TKGeomBase TKG3d)
target_include_directories(GMSHPlugin PUBLIC
  $<BUILD_INTERFACE:${GMSHPlugin_SRC_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/GMSHPlugin>
)
target_include_directories(GMSHPlugin PUBLIC ${Gmsh_PRIVATE_INCLUDE_DIRS})
target_include_directories(GMSHPlugin PRIVATE ${OpenCASCADE_INCLUDE_DIR})
if(WIN32)
    target_compile_definitions(GMSHPlugin PRIVATE GMSHPLUGIN_EXPORTS)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} GMSHPlugin)
file(GLOB GMSHPlugin_INCLUDES ${GMSHPlugin_SRC_DIR}/*.h*)
install(FILES ${GMSHPlugin_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/GMSHPlugin)

install(TARGETS ${SMESH_LIBRARIES} EXPORT SMESH-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
