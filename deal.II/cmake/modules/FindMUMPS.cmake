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

#
# TODO: mumps might link to scotch and or metis as well. Ignore this for
#       now. :-]
#

FIND_PATH(MUMPS_INCLUDE_DIRS dmumps_c.h
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

SET(MUMPS_LIBRARIES
  ${DMUMPS_LIBRARY}
  ${MUMPS_COMMON_LIBRARY}
  ${LAPACK_LIBRARIES}
  )

#
# If we can find libport.a (or similiar), link it in as well:
#
FIND_LIBRARY(PORD_LIBRARY
  NAMES port
  HINTS
    ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )
IF(NOT PORD_LIBRARY MATCHES "-NOTFOUND")
  LIST(APPEND MUMPS_LIBRARIES
    ${PORD_LIBRARY}
    )
ENDIF()

SET(MUMPS_LINKER_FLAGS
  ${LAPACK_LINKER_FLAGS}
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(MUMPS DEFAULT_MSG
  DMUMPS_LIBRARY
  MUMPS_COMMON_LIBRARY
  SCALAPACK_FOUND
  )

IF(MUMPS_FOUND)
  MARK_AS_ADVANCED(
    DMUMPS_LIBRARY
    MUMPS_COMMON_LIBRARY
    MUMPS_INCLUDE_DIRS
    PORT_LIBRARY
  )
ENDIF()

