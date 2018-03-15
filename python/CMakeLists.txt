find_package(PythonInterp 2.7 REQUIRED)
find_package(SWIG REQUIRED)
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/pyNemosys.i)
# ${CMAKE_SOURCE_DIR}/src/Mesh/meshBase.C)
set(WRAPPER_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

configure_file(${SETUP_PY_IN} ${SETUP_PY})

add_custom_target(pyNemosys ALL
                  COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install --prefix ${CMAKE_INSTALL_PREFIX}/Nemosys/python)
add_dependencies(pyNemosys Nemosys)

#Set these so python can find the wrapping. Should look into generalizing this
set( env(LD_LIBRARY_PATH ${CMAKE_BINARY_PATH}/lib/))
set( env(PYTHONPATH ${CMAKE_CURRENT_SOURCE_DIR}))



set(source1 ${CMAKE_SOURCE_DIR}/python/pyNemosys.py)
set(source2 ${CMAKE_SOURCE_DIR}/python/pyNemosys_wrap.cpp)
set(destination ${CMAKE_INSTALL_PREFIX}/Nemosys/python/src)

# writing scripts at config time so we can add them
# as commands to be executed contingent on whether pyNemosys has been built
#   - the conditionals here allow for proper make install
file(WRITE ${CMAKE_BINARY_DIR}/rename1
"if(EXISTS ${source1})
    file(RENAME ${source1} ${destination}/pyNemosys.py)
 endif()")

file(WRITE ${CMAKE_BINARY_DIR}/rename2
"if(EXISTS ${source2})
    file(RENAME ${source2} ${destination}/pyNemosys_wrap.cpp)
 endif()")

add_custom_command(TARGET pyNemosys POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/rename1
                   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/rename2)
                  

