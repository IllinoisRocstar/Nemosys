find_package(PythonInterp 3 REQUIRED) # Needed for PYTHON_EXECUTABLE
find_package(SWIG REQUIRED)

# setuptools should add Python include directories
get_target_property(WRAPPER_INCLUDE_DIRS Nemosys INTERFACE_INCLUDE_DIRECTORIES)
set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
set(SETUP_PY_IN_T "${CMAKE_CURRENT_BINARY_DIR}/setup.py.in")

# castmode so that numpy arrays can be converted to vectors/arrays
# See https://github.com/swig/swig/issues/888
set(SWIG_FLAGS "-O;-py3;-c++;-castmode")
get_target_property(SWIG_DEFS Nemosys COMPILE_DEFINITIONS)
set(SWIG_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/python)

# Could try to figure out which compiler "python setup.py build_ext" will use,
# but MSVC ignores this flag anyways
set(WRAPPER_COMPILE_FLAGS "-std=c++11")
get_target_property(WRAPPER_DEFS Nemosys INTERFACE_COMPILE_DEFINITIONS)

# Copy because python setup.py build_ext does not expose option to change where
# the resulting wrapper files (.cpp but in particular pyNemosys.py) are located
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pyNemosys.i ${CMAKE_CURRENT_BINARY_DIR}/pyNemosys.i COPYONLY)
set(WRAPPER_SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/pyNemosys.i)

# Finding vtkCommonCore library to pass to Python for linking (8.2.0)
if(VTK_ENABLE_KITS)
    NEM_get_imported_location(VTK_COMMON_CORE_LIB vtkCommon VTK_COMMON_CORE_LIB_PATH)
else()
    NEM_get_imported_location(VTK_COMMON_CORE_LIB vtkCommonCore VTK_COMMON_CORE_LIB_PATH)
endif()
if(VTK_COMMON_CORE_LIB)
    get_filename_component(VTK_COMMON_CORE_LIB ${VTK_COMMON_CORE_LIB} NAME)
    if(WIN32) # Remove the .lib suffix
        string(REGEX REPLACE "\\.lib$" "" VTK_COMMON_CORE_LIB ${VTK_COMMON_CORE_LIB})
    elseif (APPLE)
        string(REGEX REPLACE "\\.dylib$" "" VTK_COMMON_CORE_LIB ${VTK_COMMON_CORE_LIB})
        string(REGEX REPLACE "lib" "" VTK_COMMON_CORE_LIB ${VTK_COMMON_CORE_LIB})
    endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX) # Add the : in front to force exact match by ld
    set(VTK_COMMON_CORE_LIB ":${VTK_COMMON_CORE_LIB}")
endif()
message(STATUS "vtkCommonCore library ${VTK_COMMON_CORE_LIB}")

configure_file(${SETUP_PY_IN} ${SETUP_PY_IN_T})
file(GENERATE
        OUTPUT ${SETUP_PY}
        INPUT ${SETUP_PY_IN_T})

add_custom_target(pyNemosys ALL
        COMMAND "${PYTHON_EXECUTABLE}" ${SETUP_PY} build_ext --swig ${SWIG_EXECUTABLE} build)
add_dependencies(pyNemosys Nemosys)

install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" -m pip install . --user
        WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")")
