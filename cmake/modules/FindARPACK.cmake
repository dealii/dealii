## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Try to find the ARPACK library
#
# This module exports
#
#   ARPACK_LIBRARIES
#   ARPACK_LINKER_FLAGS
#   ARPACK_WITH_PARPACK
#

SET(ARPACK_DIR "" CACHE PATH "An optional hint to an ARPACK installation")
SET_IF_EMPTY(ARPACK_DIR "$ENV{ARPACK_DIR}")

SET(PARPACK_DIR "" CACHE PATH "An optional hint to a PARPACK installation")
SET_IF_EMPTY(PARPACK_DIR "$ENV{PARPACK_DIR}")

DEAL_II_FIND_LIBRARY(ARPACK_LIBRARY
  NAMES arpack
  HINTS ${ARPACK_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

IF(DEAL_II_WITH_MPI)
  GET_FILENAME_COMPONENT(_path "${ARPACK_LIBRARY}" PATH)
  DEAL_II_FIND_LIBRARY(PARPACK_LIBRARY
    NAMES parpack
    HINTS ${_path} ${ARPACK_DIR} ${PARPACK_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
    )
ELSE()
  SET(PARPACK_LIBRARY "PARPACK_LIBRARY-NOTFOUND")
ENDIF()

IF(NOT DEAL_II_ARPACK_WITH_PARPACK)
  #
  # We have to avoid an unfortunate symbol clash with "libscalapack.so" -
  # arpack happened to blindly copy a symbol name...
  #   https://github.com/opencollab/arpack-ng/issues/18
  #   https://github.com/opencollab/arpack-ng/pull/21
  #
  # Just disable parpack support if scalapack is present in Trilinos' or
  # PETSc's link interface. This can be overridden by manually setting
  # DEAL_II_ARPACK_WITH_PARPACK to true.
  #
  FOREACH(_libraries ${TRILINOS_LIBRARIES} ${PETSC_LIBRARIES})
    IF("${_libraries}" MATCHES "scalapack")
      SET(PARPACK_LIBRARY "PARPACK_LIBRARY-NOTFOUND")
    ENDIF()
  ENDFOREACH()
ENDIF()


IF(NOT PARPACK_LIBRARY MATCHES "-NOTFOUND")
  SET(ARPACK_WITH_PARPACK TRUE)
ELSE()
  SET(ARPACK_WITH_PARPACK FALSE)
ENDIF()

DEAL_II_PACKAGE_HANDLE(ARPACK
  LIBRARIES
    OPTIONAL PARPACK_LIBRARY
    REQUIRED ARPACK_LIBRARY LAPACK_LIBRARIES
    OPTIONAL MPI_C_LIBRARIES
  LINKER_FLAGS OPTIONAL LAPACK_LINKER_FLAGS
  CLEAR ARPACK_LIBRARY PARPACK_LIBRARY
  )
