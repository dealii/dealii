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
# Try to find the MUMPS library
#
# This module exports
#
#   MUMPS_INCLUDE_DIRS
#   MUMPS_LIBRARIES
#   MUMPS_LINKER_FLAGS
#

SET_IF_EMPTY(MUMPS_DIR "$ENV{MUMPS_DIR}")

INCLUDE(FindPackageHandleStandardArgs)

#
# Search for all known dependencies of MUMPS:
# (We'll rely on the user of FindMUMPS, setting up mpi *cough*)
#
FIND_PACKAGE(SCALAPACK) # which will also include lapack and blas
FIND_PACKAGE(METIS)


FIND_PATH(MUMPS_INCLUDE_DIR dmumps_c.h
  HINTS
    ${MUMPS_DIR}
  PATH_SUFFIXES
    mumps include/mumps include
  )

FIND_LIBRARY(DMUMPS_LIBRARY
  NAMES dmumps
  HINTS
    ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

FIND_LIBRARY(MUMPS_COMMON_LIBRARY
  NAMES mumps_common
  HINTS
    ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# If we can find libport.so (or similiar), link it in as well:
#
FIND_LIBRARY(PORD_LIBRARY
  NAMES pord
  HINTS
    ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )
MARK_AS_ADVANCED(PORD_LIBRARY)
IF(PORD_LIBRARY MATCHES "-NOTFOUND")
  SET(PORD_LIBRARY "")
  UNSET(PORD_LIBRARY CACHE)
ENDIF()

SET(_output ${DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${PORD_LIBRARY})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MUMPS DEFAULT_MSG
  _output # Cosmetic: Gives nice output
  DMUMPS_LIBRARY
  MUMPS_COMMON_LIBRARY
  MUMPS_INCLUDE_DIR
  SCALAPACK_FOUND
  )

MARK_AS_ADVANCED(
  DMUMPS_LIBRARY
  MUMPS_COMMON_LIBRARY
  MUMPS_INCLUDE_DIR
  PORT_LIBRARY
  )

IF(MUMPS_FOUND)
  SET(MUMPS_INCLUDE_DIRS
    ${MUMPS_INCLUDE_DIR}
    )
  SET(MUMPS_LIBRARIES
    ${DMUMPS_LIBRARY}
    ${MUMPS_COMMON_LIBRARY}
    ${PORD_LIBRARY}
    ${SCALAPACK_LIBRARIES}
    ${METIS_LIBRARIES}       # for good measure
    ${MPI_Fortran_LIBRARIES} # for good measure
    )
  SET(MUMPS_LINKER_FLAGS
    ${SCALAPACK_LINKER_FLAGS}
    )

  MARK_AS_ADVANCED(MUMPS_DIR)
ELSE()
  SET(MUMPS_DIR "" CACHE PATH
    "An optional hint to a mumps directory"
    )
ENDIF()

