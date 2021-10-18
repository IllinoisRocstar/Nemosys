# --------------------------------------------------------------------------- #
# SMESH
# --------------------------------------------------------------------------- #
set(SMESH_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Driver
file(GLOB Driver_SRC ${SMESH_SRC_DIR}/Driver/*.cxx)
add_library(Driver ${Driver_SRC})
target_link_libraries(Driver PUBLIC SMESHUtils)
target_link_libraries(Driver PRIVATE Kernel)
target_include_directories(Driver PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/Driver>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/Driver>
)
if(WIN32)
    target_compile_definitions(Driver PRIVATE MESHDRIVER_EXPORTS)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} Driver)
file(GLOB Driver_INCLUDES ${SMESH_SRC_DIR}/Driver/*.h*)
install(FILES ${Driver_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/Driver)

# DriverDAT
file(GLOB DriverDAT_SRC ${SMESH_SRC_DIR}/DriverDAT/*.cxx)
add_library(DriverDAT ${DriverDAT_SRC})
target_link_libraries(DriverDAT PUBLIC Driver)
target_link_libraries(DriverDAT PRIVATE Kernel SMDS)
target_include_directories(DriverDAT PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/DriverDAT>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverDAT>
)
if(WIN32)
    target_compile_definitions(DriverDAT PRIVATE MESHDRIVERDAT_EXPORTS)
    target_compile_definitions(DriverDAT PRIVATE UNICODE _UNICODE)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} DriverDAT)
file(GLOB DriverDAT_INCLUDES ${SMESH_SRC_DIR}/DriverDAT/*.h*)
install(FILES ${DriverDAT_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverDAT)

# DriverGMF
file(GLOB DriverGMF_SRC ${SMESH_SRC_DIR}/DriverGMF/*.cxx ${SMESH_SRC_DIR}/DriverGMF/libmesh5.c)
add_library(DriverGMF ${DriverGMF_SRC})
target_link_libraries(DriverGMF PUBLIC Driver SMDS SMESHUtils)
target_link_libraries(DriverGMF PRIVATE Kernel SMESHDS Boost::filesystem)
target_include_directories(DriverGMF PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/DriverGMF>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverGMF>
)
if(WIN32)
    target_compile_definitions(DriverGMF PRIVATE MESHDriverGMF_EXPORTS)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} DriverGMF)
file(GLOB DriverGMF_INCLUDES ${SMESH_SRC_DIR}/DriverGMF/*.h*)
install(FILES ${DriverGMF_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverGMF)

# DriverSTL
file(GLOB DriverSTL_SRC ${SMESH_SRC_DIR}/DriverSTL/*.cxx)
add_library(DriverSTL ${DriverSTL_SRC})
target_link_libraries(DriverSTL PUBLIC Driver SMDS SMESHUtils)
target_link_libraries(DriverSTL PRIVATE Kernel SMESHUtils TKernel TKSTL TKTopAlgo TKMesh)
target_include_directories(DriverSTL PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/DriverSTL>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverSTL>
)
target_include_directories(DriverSTL PRIVATE ${OpenCASCADE_INCLUDE_DIR})
if(WIN32)
    target_compile_definitions(DriverSTL PRIVATE MESHDRIVERSTL_EXPORTS)
    target_compile_definitions(DriverSTL PRIVATE UNICODE _UNICODE)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} DriverSTL)
file(GLOB DriverSTL_INCLUDES ${SMESH_SRC_DIR}/DriverSTL/*.h*)
install(FILES ${DriverSTL_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverSTL)

# DriverUNV
file(GLOB DriverUNV_SRC ${SMESH_SRC_DIR}/DriverUNV/*.cxx)
add_library(DriverUNV ${DriverUNV_SRC})
target_link_libraries(DriverUNV PUBLIC Driver SMESHDS)
target_link_libraries(DriverUNV PRIVATE Kernel SMDS)
target_include_directories(DriverUNV PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/DriverUNV>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverUNV>
)
if(WIN32)
    target_compile_definitions(DriverUNV PRIVATE MESHDRIVERUNV_EXPORTS)
    target_compile_definitions(DriverUNV PRIVATE UNICODE _UNICODE)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} DriverUNV)
file(GLOB DriverUNV_INCLUDES ${SMESH_SRC_DIR}/DriverUNV/*.h*)
install(FILES ${DriverUNV_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/DriverUNV)

# MEFISTO2
set(MEFISTO2_SRC ${SMESH_SRC_DIR}/MEFISTO2/aptrte.cxx)
include(CheckLanguage)
check_language(Fortran)
set(MEFISTO2_USING_FORTRAN OFF)
if(CMAKE_Fortran_COMPILER)
    set(MEFISTO2_USING_FORTRAN ON)
    enable_language(Fortran)
    list(APPEND MEFISTO2_SRC ${SMESH_SRC_DIR}/MEFISTO2/trte.f)
else()
    list(APPEND MEFISTO2_SRC ${SMESH_SRC_DIR}/MEFISTO2/trte.c)
endif()
add_library(MEFISTO2 ${MEFISTO2_SRC})
target_link_libraries(MEFISTO2 PUBLIC TKMath)
target_link_libraries(MEFISTO2 PRIVATE Kernel)
target_include_directories(MEFISTO2 PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/MEFISTO2>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/MEFISTO2>
  ${OpenCASCADE_INCLUDE_DIR}
)
if (WIN32)
    target_compile_definitions(MEFISTO2 PRIVATE MEFISTO2D_EXPORTS)
endif()
if(MEFISTO2_USING_FORTRAN AND CMAKE_SIZEOF_VOID_P EQUAL 8)
    target_compile_definitions(MEFISTO2 PUBLIC PCLINUX64)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} MEFISTO2)
file(GLOB MEFISTO2_INCLUDES ${SMESH_SRC_DIR}/MEFISTO2/*.h*)
install(FILES ${MEFISTO2_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/MEFISTO2)

# Controls
file(GLOB Controls_SRC ${SMESH_SRC_DIR}/Controls/*.cxx)
add_library(Controls ${Controls_SRC})
target_link_libraries(Controls PUBLIC SMESHDS SMESHUtils TKernel TKBRep TKG3d TKTopAlgo TKGeomBase TKGeomAlgo)
target_link_libraries(Controls PRIVATE Kernel Geom SMDS)
if(VTK_VERSION VERSION_GREATER_EQUAL 8.90)
  target_link_libraries(Controls PRIVATE VTK::FiltersVerdict)
else()
  target_link_libraries(Controls PRIVATE vtkFiltersVerdict)
  target_include_directories(Controls PRIVATE ${VTK_INCLUDE_DIRS})
endif()
target_include_directories(Controls PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/Controls>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/Controls>
  ${OpenCASCADE_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
)
if(WIN32)
    target_compile_definitions(Controls PRIVATE SMESHCONTROLS_EXPORTS)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} Controls)
file(GLOB Controls_INCLUDES ${SMESH_SRC_DIR}/Controls/*.h*)
install(FILES ${Controls_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/Controls)

# SMDS
file(GLOB SMDS_SRC ${SMESH_SRC_DIR}/SMDS/*.cxx)
add_library(SMDS ${SMDS_SRC})
target_link_libraries(SMDS PUBLIC Kernel)
if(VTK_VERSION VERSION_GREATER_EQUAL 8.90)
  target_link_libraries(SMDS PUBLIC VTK::CommonCore VTK::CommonDataModel)
else()
  target_link_libraries(SMDS PUBLIC vtkCommonCore vtkCommonDataModel)
  target_include_directories(SMDS PUBLIC ${VTK_INCLUDE_DIRS})
endif()
target_include_directories(SMDS PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/SMDS>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMDS>
  ${Boost_INCLUDE_DIRS}
)
if(WIN32)
    target_compile_definitions(SMDS PRIVATE SMDS_EXPORTS)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} SMDS)
file(GLOB SMDS_INCLUDES ${SMESH_SRC_DIR}/SMDS/*.h*)
install(FILES ${SMDS_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMDS)

# SMESH
file(GLOB SMESH_SRC ${SMESH_SRC_DIR}/SMESH/*.cxx)
add_library(SMESH ${SMESH_SRC})
target_link_libraries(SMESH PUBLIC Kernel SMDS SMESHDS SMESHUtils TKBRep TKMath TKG3d TKShHealing TKernel)
#target_link_libraries(SMESH PRIVATE Geom Controls Driver DriverDAT DriverSTL DriverUNV DriverGMF TKPrim TKG2d TKCDF TKMeshVS Boost::filesystem Boost::boost Boost::thread)
target_link_libraries(SMESH PRIVATE Geom Controls Driver DriverDAT DriverSTL DriverUNV DriverGMF TKPrim TKG2d TKCDF Boost::filesystem Boost::boost Boost::thread)
target_include_directories(SMESH PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/SMESH>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMESH>
  ${OpenCASCADE_INCLUDE_DIR}
)
if(WIN32)
    target_compile_definitions(SMESH PRIVATE SMESHimpl_EXPORTS)
    target_compile_definitions(SMESH PRIVATE UNICODE _UNICODE)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} SMESH)
file(GLOB SMESH_INCLUDES ${SMESH_SRC_DIR}/SMESH/*.h*)
install(FILES ${SMESH_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMESH)

# SMESHDS
file(GLOB SMESHDS_SRC ${SMESH_SRC_DIR}/SMESHDS/*.cxx)
add_library(SMESHDS ${SMESHDS_SRC})
target_link_libraries(SMESHDS PUBLIC SMDS TKBRep TKernel)
target_link_libraries(SMESHDS PRIVATE Kernel)
target_include_directories(SMESHDS PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/SMESHDS>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMESHDS>
  ${OpenCASCADE_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
)
if(WIN32)
    target_compile_definitions(SMESHDS PRIVATE SMESHDS_EXPORTS)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} SMESHDS)
file(GLOB SMESHDS_INCLUDES ${SMESH_SRC_DIR}/SMESHDS/*.h*)
install(FILES ${SMESHDS_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMESHDS)

# StdMeshers
file(GLOB StdMeshers_SRC ${SMESH_SRC_DIR}/StdMeshers/*.cxx)
add_library(StdMeshers ${StdMeshers_SRC})
target_link_libraries(StdMeshers PUBLIC Kernel SMESH SMDS SMESHUtils MEFISTO2 TKernel TKMath TKBRep TKG2d TKG3d TKShHealing)
target_link_libraries(StdMeshers PRIVATE Controls Geom TKGeomBase TKGeomAlgo TKTopAlgo TKOffset Boost::serialization)
target_include_directories(StdMeshers PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/StdMeshers>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/StdMeshers>
  ${OpenCASCADE_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
)
if(WIN32)
    target_compile_definitions(StdMeshers PRIVATE STDMESHERS_EXPORTS)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} StdMeshers)
file(GLOB StdMeshers_INCLUDES ${SMESH_SRC_DIR}/StdMeshers/*.h*)
install(FILES ${StdMeshers_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/StdMeshers)

# SMESHUtils
file(GLOB SMESHUtils_SRC ${SMESH_SRC_DIR}/SMESHUtils/*.cxx)
add_library(SMESHUtils ${SMESHUtils_SRC})
target_link_libraries(SMESHUtils PUBLIC Kernel SMDS TKBRep TKMath TKMesh)
target_link_libraries(SMESHUtils PRIVATE TKShHealing TKPrim TKernel TKBRep TKG2d TKG3d TKGeomBase TKGeomAlgo TKTopAlgo Boost::filesystem)
target_include_directories(SMESHUtils PUBLIC
  $<BUILD_INTERFACE:${SMESH_SRC_DIR}/SMESHUtils>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMESHUtils>
  ${OpenCASCADE_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
)
if(WIN32)
    target_compile_definitions(SMESHUtils PRIVATE SMESHUtils_EXPORTS)
    target_compile_definitions(SMESHUtils PRIVATE UNICODE _UNICODE)
endif()
set(SMESH_LIBRARIES ${SMESH_LIBRARIES} SMESHUtils)
file(GLOB SMESHUtils_INCLUDES ${SMESH_SRC_DIR}/SMESHUtils/*.h*)
install(FILES ${SMESHUtils_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/SMESHUtils)

install(TARGETS ${SMESH_LIBRARIES} EXPORT SMESH-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
