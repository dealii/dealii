CMAKE_MINIMUM_REQUIRED(VERSION 3.3.0)
INCLUDE(../setup_testsubproject.cmake)
PROJECT(testsuite CXX)

#
# Try to locate valgrind:
#

FIND_PROGRAM(VALGRIND_EXECUTABLE
  NAMES valgrind
  HINTS ${VALGRIND_DIR}
  PATH_SUFFIXES bin
  )

FIND_PATH(VALGRIND_INCLUDE_DIR valgrind/callgrind.h
  HINTS ${VALGRIND_DIR}
  PATH_SUFFIXES include
  )

IF( NOT VALGRIND_EXECUTABLE MATCHES "-NOTFOUND" AND
    NOT VALGRIND_INCLUDE_DIR MATCHES "-NOTFOUND" )
  MESSAGE(STATUS "Found valgrind at ${VALGRIND_EXECUTABLE}")
  SET(DEAL_II_WITH_VALGRIND TRUE)
  SET(_command_line
    "${VALGRIND_EXECUTABLE}"
    --tool=callgrind -q --combine-dumps=yes --instr-atstart=no
    --callgrind-out-file=callgrind.out
    )
ELSE()
  MESSAGE(STATUS "Could not find valgrind")
  SET(DEAL_II_WITH_VALGRIND FALSE)
ENDIF()

SET(performance_instrumentation_step_3_RUNARGS_PREFIX "${_command_line}")
SET(performance_instrumentation_step_22_RUNARGS_PREFIX "${_command_line}")

DEAL_II_PICKUP_TESTS()
