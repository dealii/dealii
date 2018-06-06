## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Try to find the NANOFLANN library
#
# This module exports
#
#   NANOFLANN_INCLUDE_DIRS
#

SET(NANOFLANN_DIR "" CACHE PATH "An optional hint to a NANOFLANN installation")
SET_IF_EMPTY(NANOFLANN_DIR "$ENV{NANOFLANN_DIR}")

DEAL_II_FIND_PATH(NANOFLANN_INCLUDE_DIR nanoflann.hpp
  HINTS ${NANOFLANN_DIR}
  PATH_SUFFIXES include
  )

DEAL_II_PACKAGE_HANDLE(NANOFLANN
  INCLUDE_DIRS REQUIRED NANOFLANN_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED NANOFLANN_INCLUDE_DIR
  CLEAR
    NANOFLANN_INCLUDE_DIR 
  )
