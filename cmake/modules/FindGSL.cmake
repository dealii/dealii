## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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
# Try to find the GSL library
#
# This module exports
#
#   GSL_LIBRARIES
#   GSL_INCLUDE_DIRS
#   GSL_LINKER_FLAGS
#

SET(GSL_DIR "" CACHE PATH "An optional hint to a GSL installation")
SET_IF_EMPTY(GSL_DIR "$ENV{GSL_DIR}")

DEAL_II_FIND_LIBRARY(GSL_LIBRARY
  NAMES gsl
  HINTS ${GSL_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_version.h
  HINTS ${GSL_DIR}
  PATH_SUFFIXES include
  )

DEAL_II_PACKAGE_HANDLE(GSL
  LIBRARIES
    REQUIRED GSL_LIBRARY
  INCLUDE_DIRS
    REQUIRED GSL_INCLUDE_DIR
  USER_INCLUDE_DIRS
    REQUIRED GSL_INCLUDE_DIR
  CLEAR GSL_LIBRARY GSL_INCLUDE_DIR
  )
