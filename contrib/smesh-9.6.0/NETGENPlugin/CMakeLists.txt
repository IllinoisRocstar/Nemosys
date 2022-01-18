# NETGENPlugin
set(NETGENPlugin_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/NETGENPlugin)

# Get sources and remove unsupported ones
file(GLOB NETGENPlugin_SRC ${NETGENPlugin_SRC_DIR}/*.cxx)
file(GLOB NETGENPlugin_Remesher_2D_SRC ${NETGENPlugin_SRC_DIR}/*Remesher_2D.cxx)
list(REMOVE_ITEM NETGENPlugin_SRC ${NETGENPlugin_Remesher_2D_SRC})
add_library(NETGENPlugin ${NETGENPlugin_SRC})
target_link_libraries(NETGENPlugin PUBLIC Kernel SMESH SMDS StdMeshers TKBRep)
target_link_libraries(NETGENPlugin PRIVATE Controls SMESHDS SMESHUtils nglib TKMesh TKMath TKGeomBase TKernel TKLCAF TKXSBase)
target_include_directories(NETGENPlugin PUBLIC
  $<BUILD_INTERFACE:${NETGENPlugin_SRC_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/NETGENPlugin>
  ${OpenCASCADE_INCLUDE_DIR}
  ${NETGEN_INCLUDE_DIR}
  ${NETGEN_INCLUDE_DIR}/include
)
target_compile_definitions(NETGENPlugin PRIVATE ${NETGEN_COMPILE_DEFINITIONS})
if(NETGEN_VERSION VERSION_GREATER_EQUAL 6)
    target_compile_definitions(NETGENPlugin PRIVATE NETGEN_V6)
    # NETGENPlugin tries to include exception.hpp in NETGENPlugin_NETGEN_3D.cxx
    # Developers of netgen probably don't intend for this to be included directly, but this approach
    # minimizes changes to source code.
    target_include_directories(NETGENPlugin PRIVATE ${NETGEN_INCLUDE_DIR}/core)
    if(NETGEN_VERSION VERSION_GREATER_EQUAL 6.2)
        target_compile_features(NETGENPlugin PRIVATE cxx_std_17)
    endif()
elseif(NETGEN_VERSION VERSION_GREATER_EQUAL 5)
    target_compile_definitions(NETGENPlugin PRIVATE NETGEN_V5)
endif()
if(WIN32)
    set_target_properties(NETGENPlugin PROPERTIES COMPILE_FLAGS "-DNETGENPLUGIN_EXPORTS")
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} NETGENPlugin)
file(GLOB NETGENPlugin_INCLUDES ${NETGENPlugin_SRC_DIR}/*.h*)
install(FILES ${NETGENPlugin_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/NETGENPlugin)

install(TARGETS ${SMESH_LIBRARIES} EXPORT SMESH-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
