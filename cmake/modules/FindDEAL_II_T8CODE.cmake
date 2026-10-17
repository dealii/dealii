## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Try to find the T8CODE library
#
# This module exports:
#   T8CODE_DIR
#   T8CODE_VERSION
#   T8CODE_VERSION_MAJOR
#   T8CODE_VERSION_MINOR
#   T8CODE_VERSION_PATCH
#   T8CODE_WITH_MPI
#

set(T8CODE_DIR "" CACHE PATH
  "An optional hint to a t8code installation/directory"
  )
set_if_empty(T8CODE_DIR "$ENV{T8CODE_DIR}")


#
# Find cmake configurations of p4est and sc.
#
# t8code requires that cmake configuration files of p4est and sc are present.
# Since FindDEAL_II_P4EST.cmake only looks for libraries and headers, we need
# to find their cmake configuration through find_package calls.
# 
# This step can be skipped once we find p4est through the cmake interface in
# FindDEAL_II_P4EST.cmake.
#
find_package(SC CONFIG
             HINTS ${P4EST_DIR}/cmake ${P4EST_DIR}/install/cmake ${T8CODE_DIR}/cmake ${T8CODE_DIR}/install/cmake)

find_package(P4EST CONFIG
             HINTS ${P4EST_DIR}/cmake ${P4EST_DIR}/install/cmake ${T8CODE_DIR}/cmake ${T8CODE_DIR}/install/cmake)

#
# Find t8code
#
find_package(T8CODE CONFIG
             HINTS ${T8CODE_DIR}/cmake ${T8CODE_DIR}/install/cmake)

if(T8CODE_FOUND)
  #
  # Adopt variable name to be consistent within deal.II.
  #
  set(T8CODE_WITH_MPI "${T8CODE_ENABLE_MPI}")

  set(_targets T8CODE::T8)
endif()

process_feature(T8CODE
  TARGETS
    REQUIRED _targets
  CLEAR
    P4EST_DIR SC_DIR
  )
