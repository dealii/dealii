find_package(deal.II 9.2.0 REQUIRED HINTS ${DEAL_II_DIR})

set(CMAKE_BUILD_TYPE ${DEAL_II_BUILD_TYPE} CACHE STRING "" FORCE)
deal_ii_initialize_cached_variables()

foreach(_var
    DIFF_DIR NUMDIFF_DIR TEST_PICKUP_REGEX TEST_TIME_LIMIT
    TEST_MPI_RANK_LIMIT TEST_THREAD_LIMIT ENABLE_PERFORMANCE_TESTS
    TESTING_ENVIRONMENT
    )
  set_if_empty(${_var} "$ENV{${_var}}")
  set(${_var} "${${_var}}" CACHE STRING "" FORCE)
endforeach()

project(TESTSUITE CXX)

#
# A custom target that does absolutely nothing. It is used in the main
# project to trigger a "make rebuild_cache" if necessary.
#
add_custom_target(regenerate)
