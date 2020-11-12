find_package(PythonInterp 2.7 REQUIRED) # Needed for PYTHON_EXECUTABLE
find_package(PythonLibs REQUIRED) # Needed for PYTHON_INCLUDE_DIRS
find_package(SWIG REQUIRED)

list(APPEND WRAPPER_INCLUDE_DIRS
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_BINARY_DIR}
        ${PYTHON_INCLUDE_DIRS})
set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY_IN_T "${CMAKE_CURRENT_BINARY_DIR}/setup.py.in")

file(COPY
        pyNemosys.i
        ${CMAKE_SOURCE_DIR}/include/meshBase.H
        ${CMAKE_SOURCE_DIR}/include/diffMesh.H
        ${CMAKE_SOURCE_DIR}/include/vtkMesh.H
        ${CMAKE_SOURCE_DIR}/include/NemDriver.H
        ${CMAKE_SOURCE_DIR}/include/TransferDriver.H
        ${CMAKE_SOURCE_DIR}/include/TransferBase.H
        ${CMAKE_SOURCE_DIR}/include/FETransfer.H
        ${CMAKE_SOURCE_DIR}/include/RefineDriver.H
        ${CMAKE_SOURCE_DIR}/include/MeshGenDriver.H
        ${CMAKE_SOURCE_DIR}/include/MeshQualityDriver.H
        ${CMAKE_SOURCE_DIR}/include/ConversionDriver.H
        ${CMAKE_SOURCE_DIR}/include/OrderOfAccuracy.H
        ${CMAKE_SOURCE_DIR}/include/meshingParams.H
        ${CMAKE_SOURCE_DIR}/include/meshGen.H
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

get_target_property(WRAPPER_COMPILE_ARGS Nemosys COMPILE_DEFINITIONS)
get_target_property(JSONCONS_COMPILE_ARGS jsoncons INTERFACE_COMPILE_DEFINITIONS)
list(APPEND WRAPPER_COMPILE_ARGS ${JSONCONS_COMPILE_ARGS})
message(STATUS "${WRAPPER_COMPILE_ARGS}")

if (ENABLE_SIMMETRIX)
    file(COPY
            pyNemosysSimmetrix.i
            ${CMAKE_SOURCE_DIR}/include/simmetrixParams.H
            ${CMAKE_SOURCE_DIR}/include/simmetrixGen.H
            ${CMAKE_SOURCE_DIR}/include/meshStitcher.H
            ${CMAKE_SOURCE_DIR}/include/meshPartitioner.H
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyNemosysSimmetrix.i)
else ()
    set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyNemosys.i)
endif ()

configure_file(${SETUP_PY_IN} ${SETUP_PY_IN_T})
file(GENERATE
        OUTPUT ${SETUP_PY}
        INPUT ${SETUP_PY_IN_T})

#configure_file(${SETUP_PY_IN} ${SETUP_PY})


add_custom_target(pyNemosys ALL
        COMMAND "${PYTHON_EXECUTABLE}" ${SETUP_PY} build)
add_dependencies(pyNemosys Nemosys)


install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" -m pip install . --user -f -O2
        WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")")
