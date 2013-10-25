## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
# Try to find the ARPACK library
#
# This module exports
#
#   ARPACK_LIBRARIES
#   ARPACK_LINKER_FLAGS
#

#
# TODO: ARPACK and mpi...
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(ARPACK_DIR "$ENV{ARPACK_DIR}")

FIND_LIBRARY(ARPACK_LIBRARY
  NAMES arpack
  HINTS
    ${ARPACK_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(ARPACK DEFAULT_MSG
  ARPACK_LIBRARY
  LAPACK_FOUND
  )

MARK_AS_ADVANCED(
  lapack_LIBRARY
  atlas_LIBRARY
  blas_LIBRARY
  ARPACK_LIBRARY
  )

IF(ARPACK_FOUND)
  SET(ARPACK_LIBRARIES
    ${ARPACK_LIBRARY}
    ${LAPACK_LIBRARIES}
    )
  SET(ARPACK_LINKER_FLAGS
    ${LAPACK_LINKER_FLAGS}
    )

  MARK_AS_ADVANCED(ARPACK_DIR)
ELSE()
  SET(ARPACK_DIR "" CACHE PATH
    "An optional hint to an ARPACK installation"
    )
ENDIF()

