## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# Try to find the GSL library
#
# This module exports
#
#   GSL_FOUND
#   GSL_LIBRARIES
#   GSL_INCLUDE_DIRS
#   GSL_LINKER_FLAGS
#   GSL_VERSION
#

#
# OK... It could be all so easy by just calling FindGSL.cmake (shipped with
# CMake around 3.2 onwards). Unfortunately this module sets up imported
# targets for the library it found (a feature). Unfortunately, portions of
# the target information seem to be cached and are thus incompatible with
# our notion of disabling and clearing a feature *sigh*.
#
# Further we support CMake from version 2.8.8 onwards and would have to do
# the manual work anyway.
#

SET(GSL_DIR "" CACHE PATH "An optional hint to a GSL installation")
SET_IF_EMPTY(GSL_DIR "$ENV{GSL_DIR}")

DEAL_II_FIND_LIBRARY(GSL_LIBRARY
  NAMES gsl
  HINTS ${GSL_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Also pick up the cblas implementation. If libgslcblas.so (or similar) is
# found we assume that gsl has to be linked against this library,
# alternatively as a fall back try known system cblas names
#
DEAL_II_FIND_LIBRARY(GSL_CBLAS_LIBRARY
  NAMES gslcblas cblas refcblas
  HINTS ${GSL_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_PATH(GSL_INCLUDE_DIR gsl/gsl_version.h
  HINTS ${GSL_DIR}
  PATH_SUFFIXES include
  )

IF(EXISTS "${GSL_INCLUDE_DIR}/gsl/gsl_version.h" )
  FILE(STRINGS "${GSL_INCLUDE_DIR}/gsl/gsl_version.h" GSL_VERSION_STRING_LINE
    REGEX "^[ \t]*#[ \t]*define[ \t]+GSL_VERSION"
    )
  STRING(REGEX REPLACE ".*([0-9].[0-9]+).*" "\\1" GSL_VERSION
    "${GSL_VERSION_STRING_LINE}"
    )
ENDIF()

DEAL_II_PACKAGE_HANDLE(GSL
  LIBRARIES
    REQUIRED GSL_LIBRARY
    OPTIONAL GSL_CBLAS_LIBRARY
  INCLUDE_DIRS
    REQUIRED GSL_INCLUDE_DIR
  USER_INCLUDE_DIRS
    REQUIRED GSL_INCLUDE_DIR
  CLEAR GSL_LIBRARY GSL_CBLAS_LIBRARY GSL_INCLUDE_DIR
  )
