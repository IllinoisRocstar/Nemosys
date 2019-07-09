cmake_minimum_required(VERSION 3.1)

message(STATUS "ACCESSING TEST DIRECTORY: ${PROJECT_SOURCE_DIR}")

# adding gtest subdir with its own cmake file
if(WIN32)
  set(BUILD_SHARED_LIBS ON)
  add_definitions(-DGTEST_HAS_TR1_TUPLE=0)
endif()
add_subdirectory(lib/gtest-1.6.0)

# include gtest library. gtest_SOURCE_DIR is availble from subdir addition
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}/include/gtest ${gtest_SOURCE_DIR})

add_custom_target(unittests)

# Set location of test scripts, test data, and run directory
file(COPY test_data
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(TEST_SCRIPTS ${CMAKE_CURRENT_SOURCE_DIR}/test_scripts)
set(TEST_DATA ${CMAKE_CURRENT_BINARY_DIR}/test_data)

# Set environment variables for testing.
if(WIN32)
  # Note: Do not indent. The indentation whitespace is added to the variable!
  set(TEST_ENV_PATH_OPTIONS
      "PATH=${CMAKE_CURRENT_BINARY_DIR};\
${CMAKE_BINARY_DIR};\
${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0;\
${CMAKE_BINARY_DIR}/contrib/madlib-1.3.0/Contrib/ANN;\
${gtest_BINARY_DIR};\
${VTK_LIB_DIR};\
${GMSH_LIB_PATH};\
${NETGEN_LIBRARY_DIR};${NETGEN_LIBRARY_DIR}/../bin;\
${ZLIB_LIB_PATH};${ZLIB_LIB_PATH}/../bin;\
${METIS_LIB_PATH};\
${CGNS_LIB_PATH};\
${HDF5_LIB_PATH};\
${EXO_LIB_PATH};${EXO_LIB_PATH}/../bin;\
${netCDF_LIB_DIR};${netCDF_LIB_DIR}/../bin;\
${SIMMETRIX_LIB_DIR};\
$ENV{PATH}")
  set(TEST_ENV_LD_OPTIONS
      "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}")
else()
  set(TEST_ENV_PATH_OPTIONS "PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}")
  # Note: Do not indent. The indentation whitespace is added to the variable!
  set(TEST_ENV_LD_OPTIONS
      "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_BINARY_DIR}:\
${VTK_LIB_DIR}:\
${NETGEN_LIBRARY_DIR}:\
${OpenCASCADE_LIBRARY_DIR}:\
$ENV{LD_LIBRARY_PATH}")
endif()

message(STATUS "Test Env Options: ${TEST_ENV_PATH_OPTIONS} ${TEST_ENV_LD_OPTIONS}")

#-------------------------------------------------------------------------------
# NEM_add_test_executable(<test_exe_name>)
#
# Adds a test executable named run<test_exe_name>Test. Requires:
#  ${TEST_SCRIPTS}/test<test_exe_name>.C  = The executable to turn into test
# Generates:
#  run<test_exe_name>Test                 = name of executable
macro(NEM_add_test_executable TEST_EXE_NAME)
  add_executable(run${TEST_EXE_NAME}Test ${TEST_SCRIPTS}/test${TEST_EXE_NAME}.C)
  target_link_libraries(run${TEST_EXE_NAME}Test gtest gtest_main Nemosys)
  add_dependencies(unittests run${TEST_EXE_NAME}Test)
endmacro()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# NEM_add_test(<test_name> <test_exe_name> <test_pwd> [arguments...])
#
# Adds a test to Nemosys named <test_name>Test using run<test_exe_name>Test
# executable. Requires:
#  run<test_exe_name>Test         = name of executable
#                                   (use NEM_add_test_executable())
#  ${TEST_DATA}/<test_pwd>        = Directory to use as working directory
#                                   to run the test
#  [arguments...]                 = Arguments passed to executable
# Generates:
#  <test_name>Test                = name of CTest test
macro(NEM_add_test TEST_NAME TEST_EXE_NAME TEST_PWD)
  add_test(NAME ${TEST_NAME}Test
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      run${TEST_EXE_NAME}Test ${ARGN}
      WORKING_DIRECTORY ${TEST_DATA}/${TEST_PWD}
  )
endmacro()
#-------------------------------------------------------------------------------

# standard tests
NEM_add_test_executable(Cubature)
NEM_add_test_executable(Conversion)
NEM_add_test_executable(OrthoPoly)
NEM_add_test_executable(PatchRecovery)
NEM_add_test_executable(Transfer)
NEM_add_test_executable(MeshGen)
NEM_add_test_executable(PNTGen)
NEM_add_test_executable(AutoVerification)
NEM_add_test_executable(Refinement)
NEM_add_test_executable(GmshMesh)

# custom built tests
if(ENABLE_EPIC)
  NEM_add_test_executable(EP16)
endif()
if(ENABLE_CFMSH)
  NEM_add_test_executable(CFM)
  NEM_add_test_executable(CFMQual)
  NEM_add_test_executable(SHM)
  NEM_add_test_executable(BlockMSH)
  NEM_add_test_executable(PackMesh)
endif()

NEM_add_test(cubatureInterp Cubature CubatureTest
    cube_refined.vtu
    cube_refinedGaussPoints.vtp
    cube_refinedGaussPointsNoData.vtp
    integrationREF.vtu
    single-hex.vtu
)

NEM_add_test(conversion Conversion ConversionTest
    case0001_ref.vtu
    case0001.msh
    hinge.vtu
    hinge.vol
    legacyVTK1.vtk
    legacyVTK2.vtk
    legacyVTK1_ref.vtu
    legacyVTK2_ref.vtu
    gorilla.vtp
    pnt-tri.pntmesh
    pnt-tri-ref.vtu
    pnt-quad.pntmesh
    pnt-quad-ref.vtu
    pnt-hex.pntmesh
    pnt-hex-ref.vtu
    pnt-mix.pntmesh
    pnt-mix-ref.vtu
)

NEM_add_test(orthoPoly OrthoPoly OrthoPolyTest F.txt)

NEM_add_test(patchRecovery PatchRecovery PatchRecoveryTest
    case0001.vtu
    testRef.vtu
    fixedWithData.vtu
)

NEM_add_test(transfer Transfer TransferTest
    pointSource.vtu
    cellSource.vtu
    target.vtu
    pntRef.vtu
    cellRef.vtu
)

NEM_add_test(meshGen MeshGen MeshGenTest
    default.json
    hingeRef.vtu
    unif.json
    hingeUnifRef.vtu
    geom.json
    hingeGeomRef.vtu
)

NEM_add_test(refinement Refinement RefinementTest
    refine_value.json
    refined_beam_value.vtu
    gold_refined_beam_value.vtu
    refine_uniform.json
    refined_beam_uniform.vtu
    gold_refined_beam_uniform.vtu
)

NEM_add_test(autoVerif AutoVerification AutoVerificationTest
    finer.vtu
    fine.vtu
    coarse.vtu
    richardson.vtu
)

NEM_add_test(pntgen PNTGen PNTGenTest
    bench1.json bench1_conv_gold.pntmesh
    bench5.json bench5_conv_gold.pntmesh
    bench6.json bench6_conv_gold.pntmesh
)

# Disable in Win for now due to CI/CD being 32 bit.
if(NOT WIN32) # TODO: Switch Win CI/CD to 64 bit to re-enable.
NEM_add_test(gmshMesh GmshMesh GmshMeshTest
    cube2.msh
    cube4.msh
    cube41.msh
)
endif()

if(ENABLE_EPIC)
  NEM_add_test(ep16 EP16 Ep16Test
      gmsh2exo.json bctest_ref.g
      ep16.json bctest_ref.dat
  )
endif()

if(ENABLE_CFMSH)
  # cfMesh tests
  NEM_add_test(cfm_small CFM cfm_sm cfmsh.json)
  NEM_add_test(cfm_2d CFM cfm_2d cfmsh.json)
  NEM_add_test(cfm_hinge CFM cfm_hinge cfmsh.json)
  NEM_add_test(cfm_bunny CFM cfm_bunny cfmsh.json)
  NEM_add_test(cfm_bunny_oct CFM cfm_bunny_oct cfmsh.json)

  # cfMesh Quality test
  NEM_add_test(cfm_qual CFMQual cfm_qual cfmsh.json cfmsh_qual.json)

  # PackMesh test
  NEM_add_test(packmesh PackMesh packmesh PackMsh_Master.json)

  # snappyHeshMesh tests
  NEM_add_test(snappy_BoxPack SHM snappy_BoxPack shmsh.json)
  NEM_add_test(snappy_Bullet SHM snappy_Bullet shmsh.json)
  NEM_add_test(snappy_Flange SHM snappy_Flange shmsh.json)

  # BlockMesh tests
  NEM_add_test(blockMsh_Box BlockMSH blockMsh_Box blockMsh_Box.json)
  NEM_add_test(blockMsh_Sphere BlockMSH blockMsh_Sphere blockMsh_Sphere.json)
  NEM_add_test(blockMsh_TCone BlockMSH blockMsh_TCone blockMsh_TCone.json)
endif()

# Python wrapper testing
if(ENABLE_PYTHON_BINDINGS)
  find_package(PythonInterp 2.7 REQUIRED)
  if(NOT PYTHONINTERP_FOUND)
    message(STATUS "Python interpreter NOT found")
  else()
    message(STATUS "Python interpreter found")
    configure_file(test_scripts/test_pyNemosys27.py ${CMAKE_CURRENT_BINARY_DIR}/../python/test_pyNemosys27.py)
    add_test(NAME testPyNemosys
        COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
            ${PYTHON_EXECUTABLE} setup.py test
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../python)
  endif()
endif()
