cmake_minimum_required(VERSION 3.13.4)
include(../scripts/setup_testsubproject.cmake)
project(testsuite CXX)

#
# Try to locate valgrind:
#

if(ENABLE_PERFORMANCE_TESTS)
  find_program(VALGRIND_EXECUTABLE
    NAMES valgrind
    HINTS ${VALGRIND_DIR}
    PATH_SUFFIXES bin
    )

  find_path(VALGRIND_INCLUDE_DIR valgrind/callgrind.h
    HINTS ${VALGRIND_DIR}
    PATH_SUFFIXES include
    )

  if( NOT VALGRIND_EXECUTABLE MATCHES "-NOTFOUND" AND
      NOT VALGRIND_INCLUDE_DIR MATCHES "-NOTFOUND" )
    message(STATUS "Found valgrind at ${VALGRIND_EXECUTABLE}")
    set(DEAL_II_WITH_VALGRIND TRUE)

    set(_command_line
      "${VALGRIND_EXECUTABLE}"
      --tool=callgrind -q --combine-dumps=yes --instr-atstart=no
      --callgrind-out-file=callgrind.out
      )
    set(performance_instrumentation_step_3_RUNARGS_PREFIX "${_command_line}")
    set(performance_instrumentation_step_22_RUNARGS_PREFIX "${_command_line}")

  else()

    message(STATUS "Could not find valgrind")
    set(DEAL_II_WITH_VALGRIND FALSE)
  endif()
endif()

deal_ii_pickup_tests()
