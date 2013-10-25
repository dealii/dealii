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
# Try to find the SCALAPACK library
#
# Used as a helper module for FindMUMPS.cmake
#
# This module exports
#
#   SCALAPACK_LIBRARIES
#   SCALAPACK_LINKER_FLAGS
#

SET_IF_EMPTY(SCALAPACK_DIR "$ENV{SCALAPACK_DIR}")
SET_IF_EMPTY(BLACS_DIR "$ENV{BLACS_DIR}")

INCLUDE(FindPackageHandleStandardArgs)

FIND_LIBRARY(SCALAPACK_LIBRARY NAMES scalapack
  HINTS
    ${SCALAPACK_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
  )

#
# Well, depending on the version of scalapack and the distribution it might
# be necessary to search for blacs, too. So we do this in a very
# probabilistic way...
#
FOREACH(_lib blacs blacsCinit blacsF77init)
  STRING(TOUPPER "${_lib}" _lib_upper)
  FIND_LIBRARY(${_lib_upper}_LIBRARY
    NAMES ${_lib} ${_lib}_MPI-LINUX-0 ${_lib}_MPI-DARWIN-0
    HINTS
      ${BLACS_DIR}
      ${SCALAPACK_DIR}
      ${SCALAPACK_DIR}/../blacs/
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib LIB
  )
  IF(NOT ${_lib_upper}_LIBRARY MATCHES "-NOTFOUND")
    LIST(APPEND BLACS_LIBRARIES
      ${${_lib_upper}_LIBRARY}
      )
  ENDIF()
ENDFOREACH()


SET(_output ${SCALAPACK_LIBRARY} ${BLACS_LIBRARIES})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCALAPACK DEFAULT_MSG
  _output # Cosmetic: Gives nice output
  SCALAPACK_LIBRARY
  LAPACK_FOUND
  )

MARK_AS_ADVANCED(
  lapack_LIBRARY
  atlas_LIBRARY
  blas_LIBRARY
  SCALAPACK_LIBRARY
  BLACS_LIBRARY
  BLACSCINIT_LIBRARY
  BLACSF77INIT_LIBRARY
  )

IF(SCALAPACK_FOUND)
  SET(SCALAPACK_LIBRARIES
    ${SCALAPACK_LIBRARY}
    ${LAPACK_LIBRARIES}
    ${BLACS_LIBRARIES}
    ${MPI_Fortran_LIBRARIES} # for good measure
    )
  SET(SCALAPACK_LINKER_FLAGS
    ${LAPACK_LINKER_FLAGS}
    )

  MARK_AS_ADVANCED(SCALAPACK_DIR BLACS_DIR)
ELSE()
  SET(SCALAPACK_DIR "" CACHE PATH
    "An optional hint to a SCALAPACK directory"
    )
  SET(BLACS_DIR "" CACHE PATH
    "An optional hint to a BLACS directory"
    )
ENDIF()

