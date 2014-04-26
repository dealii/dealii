## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# Try to find the MUMPS library
#
# This module exports
#
#   MUMPS_INCLUDE_DIRS
#   MUMPS_LIBRARIES
#   MUMPS_LINKER_FLAGS
#   MUMPS_VERSION
#   MUMPS_VERSION_MAJOR
#   MUMPS_VERSION_MINOR
#   MUMPS_VERSION_SUBMINOR
#

SET(MUMPS_DIR "" CACHE PATH "An optional hint to a mumps directory")
SET_IF_EMPTY(MUMPS_DIR "$ENV{MUMPS_DIR}")

#
# Search for all known dependencies of MUMPS:
# (We'll rely on the user of FindMUMPS, setting up mpi *cough*)
#
FIND_PACKAGE(SCALAPACK) # which will also include lapack and blas


FIND_PATH(MUMPS_INCLUDE_DIR dmumps_c.h
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES mumps include/mumps include
  )

FIND_LIBRARY(DMUMPS_LIBRARY
  NAMES dmumps
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

FIND_LIBRARY(MUMPS_COMMON_LIBRARY
  NAMES mumps_common
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# If we can find libport.so (or similiar), link it in as well:
#
FIND_LIBRARY(PORD_LIBRARY
  NAMES pord
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )
MARK_AS_ADVANCED(PORD_LIBRARY)
IF(PORD_LIBRARY MATCHES "-NOTFOUND")
  SET(PORD_LIBRARY "")
  UNSET(PORD_LIBRARY CACHE)
ENDIF()

IF(EXISTS ${MUMPS_INCLUDE_DIR}/dmumps_c.h)
  FILE(STRINGS "${MUMPS_INCLUDE_DIR}/dmumps_c.h" MUMPS_VERSION_STRING
    REGEX "#define.*MUMPS_VERSION")
  STRING(REGEX REPLACE "^.*MUMPS_VERSION.*\"(.+)\".*" "\\1"
    MUMPS_VERSION "${MUMPS_VERSION_STRING}"
    )
  STRING(REGEX REPLACE
    "([0-9]+)\\..*" "\\1" MUMPS_VERSION_MAJOR "${MUMPS_VERSION}"
    )
  STRING(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*" "\\1" MUMPS_VERSION_MINOR "${MUMPS_VERSION}"
    )
  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" MUMPS_VERSION_SUBMINOR "${MUMPS_VERSION}"
    )
ENDIF()

DEAL_II_PACKAGE_HANDLE(MUMPS
  LIBRARIES
    REQUIRED DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY
    OPTIONAL PORD_LIBRARY
    REQUIRED SCALAPACK_LIBRARIES
    OPTIONAL METIS_LIBRARIES MPI_Fortran_LIBRARIES
  INCLUDE_DIRS
    REQUIRED MUMPS_INCLUDE_DIR
  LINKER_FLAGS
    OPTIONAL SCALAPACK_LINKER_FLAGS
  )
