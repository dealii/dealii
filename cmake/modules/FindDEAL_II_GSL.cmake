## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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

set(GSL_DIR "" CACHE PATH "An optional hint to a GSL installation")
set_if_empty(GSL_DIR "$ENV{GSL_DIR}")

deal_ii_find_library(GSL_LIBRARY
  NAMES gsl
  HINTS ${GSL_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Also pick up the cblas implementation. If libgslcblas.so (or similar) is
# found we assume that gsl has to be linked against this library,
# alternatively as a fall back try known system cblas names
#
deal_ii_find_library(GSL_CBLAS_LIBRARY
  NAMES gslcblas cblas refcblas
  HINTS ${GSL_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_path(GSL_INCLUDE_DIR gsl/gsl_version.h
  HINTS ${GSL_DIR}
  PATH_SUFFIXES include
  )

if(EXISTS "${GSL_INCLUDE_DIR}/gsl/gsl_version.h" )
  file(STRINGS "${GSL_INCLUDE_DIR}/gsl/gsl_version.h" GSL_VERSION_STRING_LINE
    REGEX "^[ \t]*#[ \t]*define[ \t]+GSL_VERSION"
    )
  string(REGEX REPLACE ".*([0-9].[0-9]+).*" "\\1" GSL_VERSION
    "${GSL_VERSION_STRING_LINE}"
    )
endif()

process_feature(GSL
  LIBRARIES
    REQUIRED GSL_LIBRARY
    OPTIONAL GSL_CBLAS_LIBRARY
  INCLUDE_DIRS
    REQUIRED GSL_INCLUDE_DIR
  CLEAR GSL_LIBRARY GSL_CBLAS_LIBRARY GSL_INCLUDE_DIR
  )
