## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

if("${TRACK}" STREQUAL "")
  set(TRACK "Build Tests")
endif()
set(TEST_PICKUP_REGEX "quick_tests")
include(${CMAKE_CURRENT_LIST_DIR}/run_testsuite.cmake)
