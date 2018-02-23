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

SUBDIRLIST(TEST_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/test_data)

MESSAGE(STATUS "TEST CASES: ${TEST_DIRS}")

FIND_PROGRAM(BASH_PROGRAM bash)

SET(test_in ${CMAKE_CURRENT_SOURCE_DIR}/test_scripts/tests_in.sh)
SET(test ${CMAKE_CURRENT_SOURCE_DIR}/test_scripts/tests.sh)


CONFIGURE_FILE( ${test_in} ${test} @ONLY )

#IF (BASH_PROGRAM)
#  FOREACH(dir ${TEST_DIRS})
#    ADD_TEST(${dir} ${BASH_PROGRAM} ${test} ${dir})
#  ENDFOREACH()
#endif (BASH_PROGRAM)

SET(TEST_DATA ${CMAKE_CURRENT_LIST_DIR}/test_data)
SET(CUBATURE_TESTDIR ${TEST_DATA}/CubatureTest)
SET(CONVERSION_TESTDIR ${TEST_DATA}/ConversionTest)
SET(ORTHOPOLY_TESTDIR ${TEST_DATA}/OrthoPolyTest)
SET(PATCHRECOVERY_TESTDIR ${TEST_DATA}/PatchRecoveryTest)

ADD_TEST(NAME cubatureInterpTest COMMAND runCubatureInterpTest ${CUBATURE_TESTDIR}/cube_refined.vtu ${CUBATURE_TESTDIR}/cube_refinedGaussPoints.vtp ${CUBATURE_TESTDIR}/cube_refinedGaussPointsNoData.vtp)

ADD_TEST(NAME conversionTest COMMAND runConversionTest ${CONVERSION_TESTDIR}/case0001_ref.vtu ${CONVERSION_TESTDIR}/case0001.msh ${CONVERSION_TESTDIR}/hinge.vtu ${CONVERSION_TESTDIR}/hinge.vol)

ADD_TEST(NAME orthoPolyTest COMMAND runOrthoPolyTest ${ORTHOPOLY_TESTDIR}/F.txt)
#ADD_TEST(NAME patchRecoveryTest COMMAND runPatchRecoveryTest ${PATCHRECOVERY_TESTDIR}/cube_refined.vtu)

#Python wrapper testing

find_package(PythonInterp 3.5 REQUIRED)

if (NOT PYTHONINTERP_FOUND)
  message(STATUS "Python interpreter NOT found")
else(NOT PYTHONINTERP_FOUND)
  message(STATUS "Python interpreter found")
  ADD_TEST(NAME testPyNemosys
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_scripts/test_pyNemosys.py
           )
  set_property(TEST testPyNemosys PROPERTY ENVIRONMENT LD_LIBRARY_PATH="${CMAKE_BINARY_DIR}/lib/")
endif (NOT PYTHONINTERP_FOUND)
