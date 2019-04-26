CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(Nemosys)

MESSAGE(STATUS "ACCESSING TEST DIRECTORY: ${PROJECT_SOURCE_DIR}")

# macro to list subdirectories of a given current
MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()


# adding gtest subdir with its own cmake file
set(BUILD_SHARED_LIBS ON)
add_definitions(-DGTEST_HAS_TR1_TUPLE=0)
ADD_SUBDIRECTORY(lib/gtest-1.6.0)

# include gtest library. gtest_SOURCE_DIR is availble from subdir addition
INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR}/include/gtest ${gtest_SOURCE_DIR})

# standard tests
ADD_EXECUTABLE(runCubatureInterpTest test_scripts/testCubature.C)
ADD_EXECUTABLE(runConversionTest test_scripts/testConversion.C)
ADD_EXECUTABLE(runOrthoPolyTest test_scripts/testOrthoPoly.C)
ADD_EXECUTABLE(runPatchRecoveryTest test_scripts/testPatchRecovery.C)
ADD_EXECUTABLE(runTransferTest test_scripts/testTransfer.C)
ADD_EXECUTABLE(runMeshGenTest test_scripts/testMeshGen.C)
ADD_EXECUTABLE(runPNTGenTest test_scripts/testPNTGen.C)
ADD_EXECUTABLE(runAutoVerifTest test_scripts/testAutoVerification.C)
TARGET_LINK_LIBRARIES(runCubatureInterpTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runConversionTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runOrthoPolyTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runPatchRecoveryTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runTransferTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runMeshGenTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runPNTGenTest gtest gtest_main Nemosys)
TARGET_LINK_LIBRARIES(runAutoVerifTest gtest gtest_main Nemosys)

add_custom_target(unittests)
add_dependencies(unittests runCubatureInterpTest
                           runConversionTest
                           runOrthoPolyTest
                           runPatchRecoveryTest
                           runTransferTest
                           runMeshGenTest
                           runPNTGenTest
                           runAutoVerifTest)

# custom built test for EPIC
IF (ENABLE_EPIC)
  ADD_EXECUTABLE(runEp16Test test_scripts/testEP16.C)
  TARGET_LINK_LIBRARIES(runEp16Test gtest gtest_main Nemosys)
  add_dependencies(unittests runEp16Test)
ENDIF()

file(COPY test_data
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

SUBDIRLIST(TEST_DIRS ${CMAKE_CURRENT_BINARY_DIR}/test_data)

SET(TEST_DATA ${CMAKE_CURRENT_BINARY_DIR}/test_data)
SET(CUBATURE_TESTDIR ${TEST_DATA}/CubatureTest)
SET(CONVERSION_TESTDIR ${TEST_DATA}/ConversionTest)
SET(ORTHOPOLY_TESTDIR ${TEST_DATA}/OrthoPolyTest)
SET(PATCHRECOVERY_TESTDIR ${TEST_DATA}/PatchRecoveryTest)
SET(TRANSFER_TESTDIR ${TEST_DATA}/TransferTest)
SET(MESHGEN_TESTDIR ${TEST_DATA}/MeshGenTest)
SET(PNTGEN_TESTDIR ${TEST_DATA}/PNTGenTest)
SET(EP16_TESTDIR ${TEST_DATA}/Ep16Test)
SET(AUTOVERIF_TESTDIR ${TEST_DATA}/AutoVerificationTest)
SET(meshgendefault_in ${MESHGEN_TESTDIR}/default_in.json)
SET(meshgendefault ${MESHGEN_TESTDIR}/default.json)
SET(meshgenunif_in ${MESHGEN_TESTDIR}/unif_in.json)
SET(meshgenunif ${MESHGEN_TESTDIR}/unif.json)
SET(meshgengeom_in ${MESHGEN_TESTDIR}/geom_in.json)
SET(meshgengeom ${MESHGEN_TESTDIR}/geom.json)
SET(pntgenbench1_in ${PNTGEN_TESTDIR}/bench1_in.json)
SET(pntgenbench1 ${PNTGEN_TESTDIR}/bench1.json)
SET(pntgenbench5_in ${PNTGEN_TESTDIR}/bench5_in.json)
SET(pntgenbench5 ${PNTGEN_TESTDIR}/bench5.json)
SET(pntgenbench6_in ${PNTGEN_TESTDIR}/bench6_in.json)
SET(pntgenbench6 ${PNTGEN_TESTDIR}/bench6.json)
SET(ep16test1_in ${EP16_TESTDIR}/gmsh2exo_in.json)
SET(ep16test1 ${EP16_TESTDIR}/gmsh2exo.json)
SET(ep16test2 ${EP16_TESTDIR}/ep16.json)
CONFIGURE_FILE( ${meshgendefault_in} ${meshgendefault} @ONLY )
CONFIGURE_FILE( ${meshgenunif_in} ${meshgenunif} @ONLY )
CONFIGURE_FILE( ${meshgengeom_in} ${meshgengeom} @ONLY )
CONFIGURE_FILE( ${pntgenbench1_in} ${pntgenbench1} @ONLY )
CONFIGURE_FILE( ${pntgenbench5_in} ${pntgenbench5} @ONLY )
CONFIGURE_FILE( ${pntgenbench6_in} ${pntgenbench6} @ONLY )
CONFIGURE_FILE( ${ep16test1_in} ${ep16test1} @ONLY )

if(WIN32)
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
  set(TEST_ENV_LD_OPTIONS
"LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}:\
${VTK_LIB_DIR}:\
${NETGEN_LIB_PATH}:\
$ENV{LD_LIBRARY_PATH}")
endif()

message(STATUS "Test Env Options: ${TEST_ENV_PATH_OPTIONS} ${TEST_ENV_LD_OPTIONS}")

ADD_TEST(NAME cubatureInterpTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runCubatureInterpTest ${CUBATURE_TESTDIR}/cube_refined.vtu
                              ${CUBATURE_TESTDIR}/cube_refinedGaussPoints.vtp
                              ${CUBATURE_TESTDIR}/cube_refinedGaussPointsNoData.vtp
                              ${CUBATURE_TESTDIR}/integrationREF.vtu
                              ${CUBATURE_TESTDIR}/single-hex.vtu
)
ADD_TEST(NAME conversionTest
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
ADD_TEST(NAME orthoPolyTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runOrthoPolyTest ${ORTHOPOLY_TESTDIR}/F.txt
)
ADD_TEST(NAME patchRecoveryTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runPatchRecoveryTest ${PATCHRECOVERY_TESTDIR}/case0001.vtu
                             ${PATCHRECOVERY_TESTDIR}/testRef.vtu
                             ${PATCHRECOVERY_TESTDIR}/fixedWithData.vtu
)
ADD_TEST(NAME transferTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runTransferTest ${TRANSFER_TESTDIR}/pointSource.vtu
                        ${TRANSFER_TESTDIR}/cellSource.vtu
                        ${TRANSFER_TESTDIR}/target.vtu
                        ${TRANSFER_TESTDIR}/pntRef.vtu
                        ${TRANSFER_TESTDIR}/cellRef.vtu
)
ADD_TEST(NAME meshGenTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runMeshGenTest ${MESHGEN_TESTDIR}/default.json
                       ${MESHGEN_TESTDIR}/hingeRef.vtu
                       ${MESHGEN_TESTDIR}/unif.json
                       ${MESHGEN_TESTDIR}/hingeUnifRef.vtu
                       ${MESHGEN_TESTDIR}/geom.json
                       ${MESHGEN_TESTDIR}/hingeGeomRef.vtu
)
ADD_TEST(NAME autVerifTest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runAutoVerifTest ${AUTOVERIF_TESTDIR}/finer.vtu
                         ${AUTOVERIF_TESTDIR}/fine.vtu
                         ${AUTOVERIF_TESTDIR}/coarse.vtu
                         ${AUTOVERIF_TESTDIR}/richardson.vtu
)
ADD_TEST(NAME pntgentest
    COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
        runPNTGenTest ${PNTGEN_TESTDIR}/bench1.json ${PNTGEN_TESTDIR}/bench1_conv_gold.pntmesh
                      ${PNTGEN_TESTDIR}/bench5.json ${PNTGEN_TESTDIR}/bench5_conv_gold.pntmesh
                      ${PNTGEN_TESTDIR}/bench6.json ${PNTGEN_TESTDIR}/bench6_conv_gold.pntmesh
)
IF (ENABLE_EPIC)
  ADD_TEST(NAME ep16test
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
          runEp16Test ${ep16test1} ${EP16_TESTDIR}/bctest_ref.g
                      ${ep16test2} ${EP16_TESTDIR}/bctest_ref.dat
  )
ENDIF()

#Python wrapper testing
IF (ENABLE_PYTHON_BINDINGS)
  FIND_PACKAGE(PythonInterp 2.7 REQUIRED)
  if (NOT PYTHONINTERP_FOUND)
    message(STATUS "Python interpreter NOT found")
  else(NOT PYTHONINTERP_FOUND)
    message(STATUS "Python interpreter found")
    ADD_TEST(NAME testPyNemosys
        COMMAND  ${CMAKE_COMMAND} -E env "${TEST_ENV_PATH_OPTIONS}" "${TEST_ENV_LD_OPTIONS}"
                                         "PATH=${SWIG_DIR};$ENV{PATH}"
                                         "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR};$ENV{PYTHONPATH}"
                                         #"PYTHONPATH=${CMAKE_INSTALL_PREFIX}/Nemosys/python/src:${CMAKE_INSTALL_PREFIX}/Nemosys/python/lib/python2.7/site-packages"
            ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_scripts/test_pyNemosys27.py --topsrcdir ${CMAKE_CURRENT_SOURCE_DIR} )
    set_property(TEST testPyNemosys
                 PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib/:${CMAKE_INSTALL_PREFIX}/vtk/lib")
  endif (NOT PYTHONINTERP_FOUND)
ENDIF()
