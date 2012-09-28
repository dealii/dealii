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
# Search for the sc library, usually bundled with p4est:
#
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
  ${P4EST_LIBRARIES}
  ${SC_LIBRARIES}
  )

SET(P4EST_INCLUDE_DIRS
  ${P4EST_INCLUDE_DIR}
  ${SC_INCLUDE_DIR}
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(P4EST DEFAULT_MSG
  SC_FOUND
  P4EST_LIBRARY
  P4EST_INCLUDE_DIR
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

