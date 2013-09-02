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
# Try to find the P4EST library
#
# This module exports:
#   P4EST_LIBRARIES
#   P4EST_INCLUDE_DIRS
#   P4EST_WITH_MPI
#   P4EST_VERSION
#   P4EST_VERSION_MAJOR
#   P4EST_VERSION_MINOR
#   P4EST_VERSION_SUBMINOR
#   P4EST_VERSION_PATCH
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(P4EST_DIR "$ENV{P4EST_DIR}")
SET_IF_EMPTY(SC_DIR "$ENV{SC_DIR}")

#
# Search for the sc library, usually bundled with p4est. If no SC_DIR was
# given, take what we chose for p4est.
#

FIND_PATH(P4EST_INCLUDE_DIR p4est_config.h
  HINTS
    ${P4EST_DIR}/FAST
    ${P4EST_DIR}/DEBUG
    ${P4EST_DIR}
  PATH_SUFFIXES
    p4est include/p4est include src
  )

FIND_PATH(SC_INCLUDE_DIR sc.h
  HINTS
    ${SC_DIR}/FAST
    ${SC_DIR}/DEBUG
    ${SC_DIR}
    ${P4EST_DIR}/FAST
    ${P4EST_DIR}/DEBUG
    ${P4EST_DIR}
  PATH_SUFFIXES
    sc include/p4est include src sc/src
  )

FIND_LIBRARY(P4EST_LIBRARY_OPTIMIZED
  NAMES p4est
  HINTS
    ${P4EST_DIR}/FAST
    ${P4EST_DIR}/DEBUG
    ${P4EST_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib src
  )

FIND_LIBRARY(P4EST_LIBRARY_DEBUG
  NAMES p4est
  HINTS
    ${P4EST_DIR}/DEBUG
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib src
  )

FIND_LIBRARY(SC_LIBRARY_OPTIMIZED
  NAMES sc
  HINTS
    ${SC_DIR}/FAST
    ${SC_DIR}/DEBUG
    ${SC_DIR}
    ${P4EST_DIR}/FAST
    ${P4EST_DIR}/DEBUG
    ${P4EST_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib src sc/src
  )

FIND_LIBRARY(SC_LIBRARY_DEBUG
  NAMES sc
  HINTS
    ${SC_DIR}/DEBUG
    ${P4EST_DIR}/DEBUG
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib src sc/src
  )

SET(_output ${P4EST_LIBRARY_OPTMIZED} ${SC_LIBRARY_OPTIMIZED})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(P4EST DEFAULT_MSG
  _output # Cosmetic: Gives nice output
  P4EST_LIBRARY_OPTIMIZED
  SC_LIBRARY_OPTIMIZED
  P4EST_INCLUDE_DIR
  SC_INCLUDE_DIR
  )

MARK_AS_ADVANCED(
  P4EST_LIBRARY_OPTIMIZED
  P4EST_LIBRARY_DEBUG
  P4EST_INCLUDE_DIR
  SC_LIBRARY_OPTIMIZED
  SC_LIBRARY_DEBUG
  SC_INCLUDE_DIR
  )


IF(P4EST_FOUND)

  IF( ( "${P4EST_LIBRARY_OPTIMIZED}" STREQUAL "${P4EST_LIBRARY_DEBUG}"
        AND
        "${SC_LIBRARY_OPTIMIZED}" STREQUAL "${SC_LIBRARY_DEBUG}" )
      OR P4EST_LIBRARY_DEBUG MATCHES "-NOTFOUND"
      OR SC_LIBRARY_DEBUG MATCHES "-NOTFOUND" )
    SET(P4EST_LIBRARIES
      ${P4EST_LIBRARY_OPTIMIZED}
      ${SC_LIBRARY_OPTIMIZED}
      )
  ELSE()
    SET(P4EST_LIBRARIES
      optimized
      ${P4EST_LIBRARY_OPTIMIZED}
      ${SC_LIBRARY_OPTIMIZED}
      debug
      ${P4EST_LIBRARY_DEBUG}
      ${SC_LIBRARY_DEBUG}
      general
      )
  ENDIF()

  LIST(APPEND P4EST_LIBRARIES
    ${LAPACK_LIBRARIES} # for good measure
    ${MPI_C_LIBRARIES} # for good measure
    )
  SET(P4EST_INCLUDE_DIRS
    ${P4EST_INCLUDE_DIR}
    ${SC_INCLUDE_DIR}
    )

  #
  # Determine mpi support of p4est:
  #
  FILE(STRINGS "${P4EST_INCLUDE_DIR}/p4est_config.h" P4EST_MPI_STRING
    REGEX "#define.*P4EST_MPI 1")
  IF("${P4EST_MPI_STRING}" STREQUAL "")
    SET(P4EST_WITH_MPI FALSE)
  ELSE()

    SET(P4EST_WITH_MPI TRUE)
  ENDIF()

  #
  # Extract version numbers:
  #
  FILE(STRINGS "${P4EST_INCLUDE_DIR}/p4est_config.h" P4EST_VERSION
    REGEX "#define P4EST_VERSION \"")
  STRING(REGEX REPLACE "^.*P4EST_VERSION.*\"([0-9]+.*)\".*" "\\1"
    P4EST_VERSION "${P4EST_VERSION}"
    )
  STRING(REGEX REPLACE
    "^([0-9]+).*$" "\\1"
    P4EST_VERSION_MAJOR "${P4EST_VERSION}")
  STRING(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*$" "\\1"
    P4EST_VERSION_MINOR "${P4EST_VERSION}")
  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
    P4EST_VERSION_SUBMINOR "${P4EST_VERSION}")

  # Now for the patch number such as in 0.3.4.1. If there
  # is no patch number, then the REGEX REPLACE will fail,
  # setting P4EST_VERSION_PATCH to P4EST_VERSION. If that
  # is the case, then set the patch number to zero
  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.[0-9]+\\.([0-9]+)?.*$" "\\1"
    P4EST_VERSION_PATCH "${P4EST_VERSION}")
  IF(${P4EST_VERSION_PATCH} STREQUAL "${P4EST_VERSION}")
    SET(P4EST_VERSION_PATCH "0")
  ENDIF()


  MARK_AS_ADVANCED(P4EST_DIR)
ELSE()
  SET(P4EST_DIR "" CACHE PATH
    "An optional hint to a p4est installation/directory"
    )
ENDIF()

