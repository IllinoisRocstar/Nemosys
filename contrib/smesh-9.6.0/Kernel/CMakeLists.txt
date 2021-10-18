set(Kernel_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(Kernel_SRCS
    # Basics
    ${Kernel_SRC_DIR}/Basics/Basics_Utils.cxx
    ${Kernel_SRC_DIR}/Basics/BasicsGenericDestructor.cxx
    # Local Trace
    ${Kernel_SRC_DIR}/SALOMELocalTrace/BaseTraceCollector.cxx
    ${Kernel_SRC_DIR}/SALOMELocalTrace/FileTraceCollector.cxx
    ${Kernel_SRC_DIR}/SALOMELocalTrace/LocalTraceBufferPool.cxx
    ${Kernel_SRC_DIR}/SALOMELocalTrace/LocalTraceCollector.cxx
    # Utils
    ${Kernel_SRC_DIR}/Utils/duplicate.cxx
    ${Kernel_SRC_DIR}/Utils/OpUtil.cxx
    ${Kernel_SRC_DIR}/Utils/Utils_SALOME_Exception.cxx
    ${Kernel_SRC_DIR}/Utils/Utils_ExceptHandlers.cxx)

add_library(Kernel ${Kernel_SRCS})
target_include_directories(Kernel PUBLIC ${Boost_INCLUDE_DIRS})
target_compile_definitions(Kernel PUBLIC BOOST_DATE_TIME_NO_LIB)

target_include_directories(Kernel PUBLIC
    $<BUILD_INTERFACE:${Kernel_SRC_DIR}/Basics>
    $<BUILD_INTERFACE:${Kernel_SRC_DIR}/SALOMELocalTrace>
    $<BUILD_INTERFACE:${Kernel_SRC_DIR}/Utils>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/smesh/Kernel>
)

if(WIN32)
    target_compile_definitions(Kernel PRIVATE BASICS_EXPORTS SALOMELOCALTRACE_EXPORTS UTILS_EXPORTS)
    target_compile_definitions(Kernel PRIVATE UNICODE _UNICODE)
endif()
target_compile_definitions(Kernel PRIVATE SMESH_ONLY)

file(GLOB Kernel_INCLUDES
     ${Kernel_SRC_DIR}/Basics/*.h*
     ${Kernel_SRC_DIR}/SALOMELocalTrace/*.h*)
list(APPEND Kernel_INCLUDES
     ${Kernel_SRC_DIR}/Utils/SALOME_Utils.hxx
     ${Kernel_SRC_DIR}/Utils/OpUtil.hxx
     ${Kernel_SRC_DIR}/Utils/Utils_SALOME_Exception.hxx
     ${Kernel_SRC_DIR}/Utils/Utils_ExceptHandlers.hxx)

install(FILES ${Kernel_INCLUDES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/smesh/Kernel)

set(Kernel_LIBRARIES Kernel)

install(TARGETS ${Kernel_LIBRARIES} EXPORT SMESH-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
