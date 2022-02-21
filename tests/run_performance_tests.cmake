# A small macro
MACRO(SET_IF_EMPTY _variable)
  IF("${${_variable}}" STREQUAL "")
    SET(${_variable} ${ARGN})
  ENDIF()
ENDMACRO()

# By default we simply skip test submission to CDash
SET(TRACK "Experimental")
SET_IF_EMPTY(SKIP_SUBMISSION TRUE)

SET(CMAKE_BUILD_TYPE Release)
SET(DEAL_II_COMPILE_EXAMPLES FALSE)

SET(ENABLE_PERFORMANCE_TESTS TRUE)

#
# Determine appropriate resource limits for performance tests:
#

IF("${TESTING_ENVIRONMENT}" STREQUAL "whistler-node")
  # Managed by Matthias
  SET(CTEST_SITE "${TESTING_ENVIRONMENT}")
  SET(TESTING_ENVIRONMENT "heavy")
ENDIF()

#
# Determine appropriate resource limits for performance tests:
#

SET_IF_EMPTY(TESTING_ENVIRONMENT "light")

IF("${TESTING_ENVIRONMENT}" STREQUAL "light")
  SET_IF_EMPTY(TEST_TIME_LIMIT 600)
  SET_IF_EMPTY(TEST_MPI_RANK_LIMIT 2)
  SET_IF_EMPTY(TEST_THREAD_LIMIT 2)
ELSEIF("${TESTING_ENVIRONMENT}" STREQUAL "medium")
  SET_IF_EMPTY(TEST_TIME_LIMIT 600)
  SET_IF_EMPTY(TEST_MPI_RANK_LIMIT 8)
  SET_IF_EMPTY(TEST_THREAD_LIMIT 8)
ELSEIF("${TESTING_ENVIRONMENT}" STREQUAL "heavy")
  SET_IF_EMPTY(TEST_TIME_LIMIT 600)
  SET_IF_EMPTY(TEST_MPI_RANK_LIMIT 32)
  SET_IF_EMPTY(TEST_THREAD_LIMIT 32)
ELSE()
  MESSAGE(FATAL_ERROR
    "The variable TESTING_ENVIRONMENT was set to the invalid value "
    "»${TESTING_ENVIRONMENT}«. Valid options are light, medium, heavy.")
ENDIF()

INCLUDE(${CMAKE_CURRENT_LIST_DIR}/run_testsuite.cmake)
