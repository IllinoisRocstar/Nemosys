find_package(PythonInterp 2.7 REQUIRED)

if (PYTHONINTERP_FOUND)
    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_BINARY_DIR}/setup.py")

    set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/pyNemosys.i)
    # ${CMAKE_SOURCE_DIR}/src/Mesh/meshBase.C)
    set(WRAPPER_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    add_custom_target(pyNemosys ALL
                      COMMAND ${PYTHON_EXECUTABLE} ${SETUP_PY} install
                      )

    add_dependencies(pyNemosys Nemosys)

    #Set these so python can find the wrapping. Should look into generalizing this
    set( env(LD_LIBRARY_PATH ${CMAKE_BINARY_PATH}/lib/))
    set( env(PYTHONPATH ${CMAKE_CURRENT_SOURCE_DIR}))

endif()
