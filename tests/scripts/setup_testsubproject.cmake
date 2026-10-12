## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2013 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

find_package(deal.II 9.8.0 REQUIRED HINTS ${DEAL_II_DIR})

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
