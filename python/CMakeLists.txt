find_package(PythonInterp 2.7 REQUIRED) # Needed for PYTHON_EXECUTABLE
find_package(PythonLibs REQUIRED) # Needed for PYTHON_INCLUDE_DIRS
find_package(SWIG REQUIRED)


set(WRAPPER_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include ${PYTHON_INCLUDE_DIRS})
set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
if (ENABLE_SYMMX)
  file(COPY pyNemosysSymmx.i
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyNemosysSymmx.i)
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup_havesymmx.py.in")
else()
  file(COPY pyNemosys.i
            ${CMAKE_SOURCE_DIR}/include/meshBase.H
            ${CMAKE_SOURCE_DIR}/include/vtkMesh.H
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyNemosys.i)
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
endif()
configure_file(${SETUP_PY_IN} ${SETUP_PY})


add_custom_target(pyNemosys ALL
                  COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build)
add_dependencies(pyNemosys Nemosys)


install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install . --user -f -O2
                              WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")
