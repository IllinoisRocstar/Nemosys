find_package(PythonInterp 2.7 REQUIRED) # Needed for PYTHON_EXECUTABLE
find_package(PythonLibs REQUIRED) # Needed for PYTHON_INCLUDE_DIRS
find_package(SWIG REQUIRED)


set(WRAPPER_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include ${PYTHON_INCLUDE_DIRS})
set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
if (ENABLE_SYMMX)
  file(COPY pyNemosysSymmx.i
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyNemosysSymmx.i)
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup_havesymmx.py.in")
else()
  file(COPY pyNemosys.i
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyNemosys.i)
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
endif()
configure_file(${SETUP_PY_IN} ${SETUP_PY})


add_custom_target(pyNemosys ALL
                  COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} build --prefix ${CMAKE_INSTALL_PREFIX}/Nemosys/python)
add_dependencies(pyNemosys Nemosys)

# set these so python can find the wrapping. Should look into generalizing this
#set( env(LD_LIBRARY_PATH ${CMAKE_BINARY_PATH}/lib/))
#set( env(PYTHONPATH ${CMAKE_CURRENT_SOURCE_DIR}))

# writing scripts at config time so we can add them
# as commands to be executed contingent on whether pyNemosys has been built
#   - the conditionals here allow for proper make install
#set(source1 ${CMAKE_SOURCE_DIR}/python/pyNemosys.py)
#set(destination ${CMAKE_INSTALL_PREFIX}/Nemosys/python/lib/python2.7/site-packages)
#file(WRITE ${CMAKE_BINARY_DIR}/rename1
#"if(EXISTS ${source1})
#    file(RENAME ${source1} ${destination}/pyNemosys.py)
# endif()")
#add_custom_command(TARGET pyNemosys POST_BUILD
#                   COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/rename1)

#file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/Nemosys/python/lib)
