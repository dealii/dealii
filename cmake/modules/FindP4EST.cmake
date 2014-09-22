## ---------------------------------------------------------------------
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

SET(P4EST_DIR "" CACHE PATH
  "An optional hint to a p4est installation/directory"
  )
SET_IF_EMPTY(P4EST_DIR "$ENV{P4EST_DIR}")
SET_IF_EMPTY(SC_DIR "$ENV{SC_DIR}")

#
# Search for the sc library, usually bundled with p4est. If no SC_DIR was
# given, take what we chose for p4est.
#

DEAL_II_FIND_PATH(SC_INCLUDE_DIR sc.h
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

DEAL_II_FIND_LIBRARY(P4EST_LIBRARY_OPTIMIZED
  NAMES p4est
  HINTS ${P4EST_DIR}/FAST ${P4EST_DIR}/DEBUG ${P4EST_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib src
  )

DEAL_II_FIND_LIBRARY(SC_LIBRARY_OPTIMIZED
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

#
# Support debug variants as well:
#

DEAL_II_FIND_LIBRARY(P4EST_LIBRARY_DEBUG
  NAMES p4est
  HINTS ${P4EST_DIR}/DEBUG
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib src
  )

DEAL_II_FIND_LIBRARY(SC_LIBRARY_DEBUG
  NAMES sc
  HINTS ${SC_DIR}/DEBUG ${P4EST_DIR}/DEBUG
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib src sc/src
  )

IF( ( "${P4EST_LIBRARY_OPTIMIZED}" STREQUAL "${P4EST_LIBRARY_DEBUG}"
      AND "${SC_LIBRARY_OPTIMIZED}" STREQUAL "${SC_LIBRARY_DEBUG}" )
    OR P4EST_LIBRARY_DEBUG MATCHES "-NOTFOUND"
    OR SC_LIBRARY_DEBUG MATCHES "-NOTFOUND" )
  SET(_libraries P4EST_LIBRARY_OPTIMIZED SC_LIBRARY_OPTIMIZED)
ELSE()
  SET(_libraries
    optimized P4EST_LIBRARY_OPTIMIZED SC_LIBRARY_OPTIMIZED
    debug P4EST_LIBRARY_DEBUG SC_LIBRARY_DEBUG
    general
    )
ENDIF()


DEAL_II_FIND_PATH(P4EST_INCLUDE_DIR p4est_config.h
  HINTS ${P4EST_DIR}/FAST ${P4EST_DIR}/DEBUG ${P4EST_DIR}
  PATH_SUFFIXES p4est include/p4est include src
  )

IF(EXISTS ${P4EST_INCLUDE_DIR}/p4est_config.h)
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
  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
    P4EST_VERSION_PATCH "${P4EST_VERSION}")

  #
  # We cannot rely on the fact that SUBMINOR or PATCH are defined.
  # Nevertheless, we need a full version number for our preprocessor macros
  # to work. If the p4est version number is only of the form x.y instead of
  # a.b.c.d, then the last two REGEX_REPLACE calls above will have failed
  # because the regular expression didn't match the version string,
  # and P4EST_VERSION_SUBMINOR and P4EST_VERSION_PATCH will either be
  # empty or be the full version string. In those cases, set those numbers
  # to 0 if necessary.
  #
  IF("${P4EST_VERSION_SUBMINOR}" MATCHES "^(|${P4EST_VERSION})$")
    SET(P4EST_VERSION_SUBMINOR "0")
  ENDIF()

  IF("${P4EST_VERSION_PATCH}" MATCHES "^(|${P4EST_VERSION})$")
    SET(P4EST_VERSION_PATCH "0")
  ENDIF()
ENDIF()

DEAL_II_PACKAGE_HANDLE(P4EST
  LIBRARIES
    REQUIRED ${_libraries}
    OPTIONAL LAPACK_LIBRARIES MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED P4EST_INCLUDE_DIR SC_INCLUDE_DIR
  USER_INCLUDE_DIRS
    REQUIRED P4EST_INCLUDE_DIR SC_INCLUDE_DIR
  CLEAR
    SC_INCLUDE_DIR P4EST_LIBRARY_OPTIMIZED SC_LIBRARY_OPTIMIZED
    P4EST_LIBRARY_DEBUG SC_LIBRARY_DEBUG P4EST_INCLUDE_DIR
  )
