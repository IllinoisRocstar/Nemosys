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

IF (BASH_PROGRAM)
  FOREACH(dir ${TEST_DIRS})
    ADD_TEST(${dir} ${BASH_PROGRAM} ${test} ${dir})
  ENDFOREACH()
endif (BASH_PROGRAM)

