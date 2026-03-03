## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
#   T8CODE_LIBRARIES
#   T8CODE_INCLUDE_DIR
#   T8CODE_WITH_MPI
#

set(T8CODE_DIR "" CACHE PATH
  "An optional hint to a t8code installation/directory"
  )
set_if_empty(T8CODE_DIR "$ENV{T8CODE_DIR}")



## add additional hints for different ways to install t8code/p4est/sc
find_package(SC CONFIG
             HINTS ${T8CODE_DIR}/cmake ${T8CODE_DIR}/install/cmake)

find_package(P4EST CONFIG
             HINTS ${T8CODE_DIR}/cmake ${T8CODE_DIR}/install/cmake)

find_package(T8CODE CONFIG
             HINTS ${T8CODE_DIR}/cmake ${T8CODE_DIR}/install/cmake)

  if(${T8CODE_ENABLE_MPI})
    message(STATUS "Found MPI")
    set(T8CODE_WITH_MPI TRUE)
  else()
    message(STATUS "NOT Found MPI")
  endif()

 
  deal_ii_find_path(T8CODE_INCLUDE_DIR t8.h
  HINTS ${T8CODE_DIR}/.. ${T8CODE_DIR}
  PATH_SUFFIXES include
  )

add_definitions(-DT8_CMAKE_BUILD)
set(_targets T8CODE::T8)

process_feature(T8CODE
  TARGETS REQUIRED _targets
  LIBRARIES OPTIONAL LAPACK_LIBRARIES MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED T8CODE_INCLUDE_DIR
  )
