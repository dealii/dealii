#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Try to find the P4EST library
#
# This module exports:
#   P4EST_LIBRARIES
#   P4EST_INCLUDE_DIRS
#   P4EST_WITH_MPI
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(P4EST_DIR "$ENV{P4EST_DIR}")

#
# We used to recommend installing p4est with a custom script that
# compiled p4est twice, once in debug and once in optimized mode.
# the installation would then have happened into directories
# $P4EST_DIR/DEBUG and $P4EST_DIR/FAST. If we can find these
# two directories, then use the FAST directory rather than trying
# to figure out how we can build deal.II against the two libraries
# depending on whether we are in debug or optimized mode.
#
IF (P4EST_DIR
    AND
    EXISTS ${P4EST_DIR}/DEBUG
    AND
    EXISTS ${P4EST_DIR}/FAST)
  MESSAGE(STATUS "Found old-style p4est directory layout")
  SET (P4EST_DIR ${P4EST_DIR}/FAST)
ENDIF()


#
# Search for the sc library, usually bundled with p4est. If no SC_DIR was
# given, take what we chose for p4est.
#
SET_IF_EMPTY(SC_DIR "$ENV{SC_DIR}")
SET_IF_EMPTY(SC_DIR "${P4EST_DIR}")
FIND_PACKAGE(SC)

FIND_PATH(P4EST_INCLUDE_DIR p4est.h
  HINTS
    ${P4EST_DIR}
  PATH_SUFFIXES
    p4est include/p4est include src
  )

FIND_LIBRARY(P4EST_LIBRARY
  NAMES p4est
  HINTS
    ${P4EST_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib src
  )

SET(P4EST_LIBRARIES
  ${P4EST_LIBRARY}
  ${SC_LIBRARY}
  )

SET(P4EST_INCLUDE_DIRS
  ${P4EST_INCLUDE_DIR}
  ${SC_INCLUDE_DIR}
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(P4EST DEFAULT_MSG
  P4EST_LIBRARY
  P4EST_INCLUDE_DIR
  SC_FOUND
  )

IF(P4EST_FOUND)
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

  MARK_AS_ADVANCED(
    P4EST_LIBRARY
    P4EST_INCLUDE_DIR
    P4EST_DIR
  )
ELSE()
  SET(P4EST_DIR "" CACHE STRING
    "An optional hint to a p4est installation/directory"
    )
ENDIF()

