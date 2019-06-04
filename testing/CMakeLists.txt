cmake_minimum_required(VERSION 3.1)

message(STATUS "ACCESSING TEST DIRECTORY: ${PROJECT_SOURCE_DIR}")

# macro to list subdirectories of a given current
macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
      list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

# adding gtest subdir with its own cmake file
if(WIN32)
  set(BUILD_SHARED_LIBS ON)
  add_definitions(-DGTEST_HAS_TR1_TUPLE=0)
endif(WIN32)
add_subdirectory(lib/gtest-1.6.0)

# include gtest library. gtest_SOURCE_DIR is availble from subdir addition
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}/include/gtest ${gtest_SOURCE_DIR})

# standard tests
add_executable(runCubatureInterpTest test_scripts/testCubature.C)
add_executable(runConversionTest test_scripts/testConversion.C)
add_executable(runOrthoPolyTest test_scripts/testOrthoPoly.C)
add_executable(runPatchRecoveryTest test_scripts/testPatchRecovery.C)
add_executable(runTransferTest test_scripts/testTransfer.C)
add_executable(runMeshGenTest test_scripts/testMeshGen.C)
add_executable(runPNTGenTest test_scripts/testPNTGen.C)
add_executable(runAutoVerifTest test_scripts/testAutoVerification.C)
add_executable(runRefinementTest test_scripts/testRefinement.C)
add_executable(runGmshMeshTest test_scripts/testGmshMesh.C)

target_link_libraries(runCubatureInterpTest gtest gtest_main Nemosys)
target_link_libraries(runConversionTest gtest gtest_main Nemosys)
target_link_libraries(runOrthoPolyTest gtest gtest_main Nemosys)
target_link_libraries(runPatchRecoveryTest gtest gtest_main Nemosys)
target_link_libraries(runTransferTest gtest gtest_main Nemosys)
target_link_libraries(runMeshGenTest gtest gtest_main Nemosys)
target_link_libraries(runPNTGenTest gtest gtest_main Nemosys)
target_link_libraries(runAutoVerifTest gtest gtest_main Nemosys)
target_link_libraries(runRefinementTest gtest gtest_main Nemosys)
target_link_libraries(runGmshMeshTest gtest gtest_main Nemosys)

add_custom_target(unittests)
add_dependencies(unittests runCubatureInterpTest
                           runConversionTest
                           runOrthoPolyTest
                           runPatchRecoveryTest
                           runTransferTest
                           runMeshGenTest
                           runPNTGenTest
                           runAutoVerifTest
                           runRefinementTest
                           runGmshMeshTest
)

# custom built tests
if(ENABLE_EPIC)
  add_executable(runEp16Test test_scripts/testEP16.C)
  target_link_libraries(runEp16Test gtest gtest_main Nemosys)
  add_dependencies(unittests runEp16Test)
endif()
if(ENABLE_CFMSH)
  add_executable(runCFMTest test_scripts/testCFM.C)
  target_link_libraries(runCFMTest gtest gtest_main Nemosys)
  add_dependencies(unittests runCFMTest)

  add_executable(runCFMQualTest test_scripts/testCFMQual.C)
  target_link_libraries(runCFMQualTest gtest gtest_main Nemosys)
  add_dependencies(unittests runCFMQualTest)
endif()

file(COPY test_data
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

SUBDIRLIST(TEST_DIRS ${CMAKE_CURRENT_BINARY_DIR}/test_data)

set(TEST_DATA ${CMAKE_CURRENT_BINARY_DIR}/test_data)
set(CUBATURE_TESTDIR ${TEST_DATA}/CubatureTest)
set(CONVERSION_TESTDIR ${TEST_DATA}/ConversionTest)
set(ORTHOPOLY_TESTDIR ${TEST_DATA}/OrthoPolyTest)
set(PATCHRECOVERY_TESTDIR ${TEST_DATA}/PatchRecoveryTest)
set(TRANSFER_TESTDIR ${TEST_DATA}/TransferTest)
set(MESHGEN_TESTDIR ${TEST_DATA}/MeshGenTest)
set(PNTGEN_TESTDIR ${TEST_DATA}/PNTGenTest)
set(EP16_TESTDIR ${TEST_DATA}/Ep16Test)
set(AUTOVERIF_TESTDIR ${TEST_DATA}/AutoVerificationTest)
set(REFINEMENT_TESTDIR ${TEST_DATA}/RefinementTest)
set(CFM_SM_TESTDIR ${TEST_DATA}/cfm_sm)
set(CFM_2D_TESTDIR ${TEST_DATA}/cfm_2d)
set(CFM_HINGE_TESTDIR ${TEST_DATA}/cfm_hinge)
set(CFM_BUNNY_TESTDIR ${TEST_DATA}/cfm_bunny)
set(CFM_BUNNY_OCT_TESTDIR ${TEST_DATA}/cfm_bunny_oct)
set(CFM_QUAL_TESTDIR ${TEST_DATA}/cfm_qual)
set(GMSHMESH_TESTDIR ${TEST_DATA}/GmshMeshTest)

set(meshgendefault_in ${MESHGEN_TESTDIR}/default_in.json)
set(meshgendefault ${MESHGEN_TESTDIR}/default.json)
set(meshgenunif_in ${MESHGEN_TESTDIR}/unif_in.json)
set(meshgenunif ${MESHGEN_TESTDIR}/unif.json)
set(meshgengeom_in ${MESHGEN_TESTDIR}/geom_in.json)
set(meshgengeom ${MESHGEN_TESTDIR}/geom.json)
set(pntgenbench1_in ${PNTGEN_TESTDIR}/bench1_in.json)
set(pntgenbench1 ${PNTGEN_TESTDIR}/bench1.json)
set(pntgenbench5_in ${PNTGEN_TESTDIR}/bench5_in.json)
set(pntgenbench5 ${PNTGEN_TESTDIR}/bench5.json)
set(pntgenbench6_in ${PNTGEN_TESTDIR}/bench6_in.json)
set(pntgenbench6 ${PNTGEN_TESTDIR}/bench6.json)
set(ep16test1_in ${EP16_TESTDIR}/gmsh2exo_in.json)
set(ep16test1 ${EP16_TESTDIR}/gmsh2exo.json)
set(ep16test2 ${EP16_TESTDIR}/ep16.json)
set(refinementValue_in ${REFINEMENT_TESTDIR}/refine_value_in.json)
set(refinementValue ${REFINEMENT_TESTDIR}/refine_value.json)
set(refinementUniform_in ${REFINEMENT_TESTDIR}/refine_uniform_in.json)
set(refinementUniform ${REFINEMENT_TESTDIR}/refine_uniform.json)
set(cfm_sm_in ${CFM_SM_TESTDIR}/cfmsh_in.json)
set(cfm_sm ${CFM_SM_TESTDIR}/cfmsh.json)
set(cfm_2d_in ${CFM_2D_TESTDIR}/cfmsh_in.json)
set(cfm_2d ${CFM_2D_TESTDIR}/cfmsh.json)
set(cfm_hinge_in ${CFM_HINGE_TESTDIR}/cfmsh_in.json)
set(cfm_hinge ${CFM_HINGE_TESTDIR}/cfmsh.json)
set(cfm_bunny_in ${CFM_BUNNY_TESTDIR}/cfmsh_in.json)
set(cfm_bunny ${CFM_BUNNY_TESTDIR}/cfmsh.json)
set(cfm_bunny_oct_in ${CFM_BUNNY_OCT_TESTDIR}/cfmsh_in.json)
set(cfm_bunny_oct ${CFM_BUNNY_OCT_TESTDIR}/cfmsh.json)
set(cfm_qual_gen_in ${CFM_QUAL_TESTDIR}/cfmsh_in.json)
set(cfm_qual_gen ${CFM_QUAL_TESTDIR}/cfmsh.json)
set(cfm_qual ${CFM_QUAL_TESTDIR}/cfmsh_qual.json)

configure_file(${meshgendefault_in} ${meshgendefault} @ONLY)
configure_file(${meshgenunif_in} ${meshgenunif} @ONLY)
configure_file(${meshgengeom_in} ${meshgengeom} @ONLY)
configure_file(${pntgenbench1_in} ${pntgenbench1} @ONLY)
configure_file(${pntgenbench5_in} ${pntgenbench5} @ONLY)
configure_file(${pntgenbench6_in} ${pntgenbench6} @ONLY)
configure_file(${ep16test1_in} ${ep16test1} @ONLY)
configure_file(${refinementValue_in} ${refinementValue} @ONLY)
configure_file(${refinementUniform_in} ${refinementUniform} @ONLY)
configure_file(${cfm_sm_in} ${cfm_sm} @ONLY)
configure_file(${cfm_2d_in} ${cfm_2d} @ONLY)
configure_file(${cfm_hinge_in} ${cfm_hinge} @ONLY)
configure_file(${cfm_bunny_in} ${cfm_bunny} @ONLY)
configure_file(${cfm_bunny_oct_in} ${cfm_bunny_oct} @ONLY)
configure_file(${cfm_qual_gen_in} ${cfm_qual_gen} @ONLY)

if(WIN32)
  # Note: Do not indent. The indentation whitespace is added to the variable!
  set(TEST_ENV_PATH_OPTIONS
"PATH=${CMAKE_BINARY_DIR};\
${gtest_BINARY_DIR};\
${MADLIB_LIB_PATH};\
${VTK_LIB_DIR};\
${GMSH_LIB_PATH};\
${NETGEN_LIB_PATH};\
${NETGEN_LIB_PATH}/../bin;\
${ZLIB_LIB_PATH};\
${ZLIB_LIB_PATH}/../bin;\
${ANN_LIB_PATH};\
${METIS_LIB_PATH};\
${CGNS_LIB_PATH};\
$ENV{PATH}")
else()
  set(TEST_ENV_PATH_OPTIONS "PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PATH}")
  # Note: Do not indent. The indentation whitespace is added to the variable!
  set(TEST_ENV_LD_OPTIONS
"LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:\
${CMAKE_CURRENT_BINARY_DIR}/..:\
${VTK_LIB_DIR}:\
${NETGEN_LIB_PATH}:\
${OpenCASCADE_LIBRARY_DIR}:\
$ENV{LD_LIBRARY_PATH}")
endif()

message(STATUS "Test Env Options: ${TEST_ENV_PATH_OPTIONS} ${TEST_ENV_LD_OPTIONS}")

add_test(NAME cubatureInterpTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runCubatureInterpTest ${CUBATURE_TESTDIR}/cube_refined.vtu
                              ${CUBATURE_TESTDIR}/cube_refinedGaussPoints.vtp
                              ${CUBATURE_TESTDIR}/cube_refinedGaussPointsNoData.vtp
                              ${CUBATURE_TESTDIR}/integrationREF.vtu
                              ${CUBATURE_TESTDIR}/single-hex.vtu
)
add_test(NAME conversionTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runConversionTest ${CONVERSION_TESTDIR}/case0001_ref.vtu
                          ${CONVERSION_TESTDIR}/case0001.msh
                          ${CONVERSION_TESTDIR}/hinge.vtu
                          ${CONVERSION_TESTDIR}/hinge.vol
                          ${CONVERSION_TESTDIR}/legacyVTK1.vtk
                          ${CONVERSION_TESTDIR}/legacyVTK2.vtk
                          ${CONVERSION_TESTDIR}/legacyVTK1_ref.vtu
                          ${CONVERSION_TESTDIR}/legacyVTK2_ref.vtu
                          ${CONVERSION_TESTDIR}/gorilla.vtp
                          ${CONVERSION_TESTDIR}/pnt-tri.pntmesh
                          ${CONVERSION_TESTDIR}/pnt-tri-ref.vtu
                          ${CONVERSION_TESTDIR}/pnt-quad.pntmesh
                          ${CONVERSION_TESTDIR}/pnt-quad-ref.vtu
                          ${CONVERSION_TESTDIR}/pnt-hex.pntmesh
                          ${CONVERSION_TESTDIR}/pnt-hex-ref.vtu
                          ${CONVERSION_TESTDIR}/pnt-mix.pntmesh
                          ${CONVERSION_TESTDIR}/pnt-mix-ref.vtu
)
add_test(NAME orthoPolyTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runOrthoPolyTest ${ORTHOPOLY_TESTDIR}/F.txt
)
add_test(NAME patchRecoveryTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runPatchRecoveryTest ${PATCHRECOVERY_TESTDIR}/case0001.vtu
                             ${PATCHRECOVERY_TESTDIR}/testRef.vtu
                             ${PATCHRECOVERY_TESTDIR}/fixedWithData.vtu
)
add_test(NAME transferTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runTransferTest ${TRANSFER_TESTDIR}/pointSource.vtu
                        ${TRANSFER_TESTDIR}/cellSource.vtu
                        ${TRANSFER_TESTDIR}/target.vtu
                        ${TRANSFER_TESTDIR}/pntRef.vtu
                        ${TRANSFER_TESTDIR}/cellRef.vtu
)
add_test(NAME meshGenTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runMeshGenTest ${MESHGEN_TESTDIR}/default.json
                       ${MESHGEN_TESTDIR}/hingeRef.vtu
                       ${MESHGEN_TESTDIR}/unif.json
                       ${MESHGEN_TESTDIR}/hingeUnifRef.vtu
                       ${MESHGEN_TESTDIR}/geom.json
                       ${MESHGEN_TESTDIR}/hingeGeomRef.vtu
)
add_test(NAME refinementTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runRefinementTest ${REFINEMENT_TESTDIR}/refine_value.json
                          ${REFINEMENT_TESTDIR}/refined_beam_value.vtu
                          ${REFINEMENT_TESTDIR}/gold_refined_beam_value.vtu
                          ${REFINEMENT_TESTDIR}/refine_uniform.json
                          ${REFINEMENT_TESTDIR}/refined_beam_uniform.vtu
                          ${REFINEMENT_TESTDIR}/gold_refined_beam_uniform.vtu
)
add_test(NAME autVerifTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runAutoVerifTest ${AUTOVERIF_TESTDIR}/finer.vtu
                         ${AUTOVERIF_TESTDIR}/fine.vtu
                         ${AUTOVERIF_TESTDIR}/coarse.vtu
                         ${AUTOVERIF_TESTDIR}/richardson.vtu
)
add_test(NAME pntgentest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runPNTGenTest ${PNTGEN_TESTDIR}/bench1.json ${PNTGEN_TESTDIR}/bench1_conv_gold.pntmesh
                      ${PNTGEN_TESTDIR}/bench5.json ${PNTGEN_TESTDIR}/bench5_conv_gold.pntmesh
                      ${PNTGEN_TESTDIR}/bench6.json ${PNTGEN_TESTDIR}/bench6_conv_gold.pntmesh
)
add_test(NAME gmshMeshTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runGmshMeshTest ${GMSHMESH_TESTDIR}/cube2.msh
                        ${GMSHMESH_TESTDIR}/cube4.msh
                        ${GMSHMESH_TESTDIR}/cube41.msh
)
if(ENABLE_EPIC)
  add_test(NAME ep16test
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
          runEp16Test ${ep16test1} ${EP16_TESTDIR}/bctest_ref.g
                      ${ep16test2} ${EP16_TESTDIR}/bctest_ref.dat
  )
endif()
if(ENABLE_CFMSH)
  add_test(NAME cfm_small
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      runCFMTest ${cfm_sm}
      WORKING_DIRECTORY ${CFM_SM_TESTDIR}
  )
  add_test(NAME cfm_2d
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      runCFMTest ${cfm_2d}
      WORKING_DIRECTORY ${CFM_2D_TESTDIR}
  )
  add_test(NAME cfm_hinge
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      runCFMTest ${cfm_hinge}
      WORKING_DIRECTORY ${CFM_HINGE_TESTDIR}
  )
  add_test(NAME cfm_bunny
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      runCFMTest ${cfm_bunny}
      WORKING_DIRECTORY ${CFM_BUNNY_TESTDIR}
  )
  add_test(NAME cfm_bunny_oct
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      runCFMTest ${cfm_bunny_oct}
      WORKING_DIRECTORY ${CFM_BUNNY_OCT_TESTDIR}
  )
  add_test(NAME cfm_qual
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
      runCFMQualTest ${cfm_qual_gen} ${cfm_qual}
      WORKING_DIRECTORY ${CFM_QUAL_TESTDIR}
  )
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
