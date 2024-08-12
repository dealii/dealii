# A small macro
macro(set_if_empty _variable)
  if("${${_variable}}" STREQUAL "")
    set(${_variable} ${ARGN})
  endif()
endmacro()

# By default we simply skip test submission to CDash
set(TRACK "Experimental")
set_if_empty(SKIP_SUBMISSION TRUE)

set(CMAKE_BUILD_TYPE Release)

set(ENABLE_PERFORMANCE_TESTS TRUE)

#
# Determine appropriate resource limits for performance tests:
#

if("${TESTING_ENVIRONMENT}" STREQUAL "whistler-node")
  # Managed by Matthias
  set(CTEST_SITE "${TESTING_ENVIRONMENT}")
  set(TESTING_ENVIRONMENT "heavy")
endif()

#
# Determine appropriate resource limits for performance tests:
#

set_if_empty(TESTING_ENVIRONMENT "light")

if("${TESTING_ENVIRONMENT}" STREQUAL "light")
  set_if_empty(TEST_TIME_LIMIT 600)
  set_if_empty(TEST_MPI_RANK_LIMIT 2)
  set_if_empty(TEST_THREAD_LIMIT 2)
elseif("${TESTING_ENVIRONMENT}" STREQUAL "medium")
  set_if_empty(TEST_TIME_LIMIT 600)
  set_if_empty(TEST_MPI_RANK_LIMIT 8)
  set_if_empty(TEST_THREAD_LIMIT 8)
elseif("${TESTING_ENVIRONMENT}" STREQUAL "heavy")
  set_if_empty(TEST_TIME_LIMIT 600)
  set_if_empty(TEST_MPI_RANK_LIMIT 32)
  set_if_empty(TEST_THREAD_LIMIT 32)
else()
  message(FATAL_ERROR
    "The variable TESTING_ENVIRONMENT was set to the invalid value "
    "»${TESTING_ENVIRONMENT}«. Valid options are light, medium, heavy.")
endif()

include(${CMAKE_CURRENT_LIST_DIR}/run_testsuite.cmake)
