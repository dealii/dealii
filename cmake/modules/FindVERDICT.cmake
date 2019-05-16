## ---------------------------------------------------------------------
##
## Copyright (C) 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find VERDICT
#
# This module exports
#
#   VERDICT_INCLUDE_DIR
#   VERDICT_LIBRARY
#

SET(VERDICT_DIR "" CACHE PATH "An optional hint to a VERDICT installation")
SET_IF_EMPTY(VERDICT_DIR "$ENV{VERDICT_DIR}")

DEAL_II_FIND_PATH(VERDICT_INCLUDE_DIR
  NAMES verdict.h
  HINTS ${VERDICT_DIR}
  PATH_SUFFIXES include
  )

DEAL_II_FIND_LIBRARY(VERDICT_LIBRARY
  NAMES verdict
  HINTS ${VERDICT_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )


DEAL_II_PACKAGE_HANDLE(VERDICT
  LIBRARIES 
    REQUIRED VERDICT_LIBRARY
  INCLUDE_DIRS REQUIRED VERDICT_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED VERDICT_INCLUDE_DIR
  CLEAR VERDICT_INCLUDE_DIR VERDICT_LIBRARY
  )
