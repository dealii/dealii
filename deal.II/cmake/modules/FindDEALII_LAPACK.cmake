## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2013 by the deal.II authors
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
# This module is a wrapper around the FindLAPACK.cmake module provided by
# CMake.
#
# This module exports
#
#   LAPACK_FOUND
#   LAPACK_LIBRARIES
#   LAPACK_LINKER_FLAGS
#   BLAS_FOUND
#   BLAS_LIBRARIES
#   BLAS_LINKER_FLAGS
#



#
# We have to use a trick with CMAKE_PREFIX_PATH to make LAPACK_DIR and
# BLAS_DIR work...
#
SET_IF_EMPTY(BLAS_DIR "$ENV{BLAS_DIR}")
SET_IF_EMPTY(LAPACK_DIR "$ENV{LAPACK_DIR}")

SET(_cmake_prefix_path_backup "${CMAKE_PREFIX_PATH}")

SET(CMAKE_PREFIX_PATH ${BLAS_DIR} ${LAPACK_DIR} ${_cmake_prefix_path_backup})

FIND_PACKAGE(BLAS)

SET(CMAKE_PREFIX_PATH ${LAPACK_DIR} ${_cmake_prefix_path_backup})

FIND_PACKAGE(LAPACK)

SET(CMAKE_PREFIX_PATH ${_cmake_prefix_path_backup})

MARK_AS_ADVANCED(
  atlas_LIBRARY atlcblas_LIBRARY atllapack_LIBRARY blas_LIBRARY
  eigen_blas_LIBRARY f77blas_LIBRARY gslcblas_LIBRARY lapack_LIBRARY
  m_LIBRARY ptf77blas_LIBRARY ptlapack_LIBRARY refblas_LIBRARY
  reflapack_LIBRARY
  )


IF(LAPACK_FOUND)
  SET(DEALII_LAPACK_FOUND TRUE)

  #
  # So, well... LAPACK_LINKER_FLAGS and LAPACK_LIBRARIES should contain the
  # complete link interface. But for invalid user overrides we include
  # BLAS_LIBRARIES and BLAS_LINKER_FLAGS as well..
  #
  IF(NOT LAPACK_LINKER_FLAGS MATCHES "${BLAS_LINKER_FLAGS}")
    MESSAGE(STATUS
      "Manually adding BLAS_LINKER_FLAGS to LAPACK_LINKER_FLAGS"
      )
    ADD_FLAGS(LAPACK_LINKER_FLAGS "${BLAS_LINKER_FLAGS}")
  ENDIF()
  IF(NOT "${LAPACK_LIBRARIES}" MATCHES "${BLAS_LIBRARIES}")
    MESSAGE(STATUS
      "Manually adding BLAS_LIBRARIES to LAPACK_LIBRARIES"
      )
    LIST(APPEND LAPACK_LIBRARIES ${BLAS_LIBRARIES})
  ENDIF()

  #
  # Well, in case of static archives we have to manually pick up the
  # complete link interface. *sigh*
  #
  # If CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is not available, do it
  # unconditionally for the most common case (gfortran).
  #
  SET(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  SET_IF_EMPTY(_fortran_libs gfortran m quadmath c)

  FOREACH(_lib ${_fortran_libs})
    FIND_SYSTEM_LIBRARY(${_lib}_LIBRARY NAMES ${_lib})
    MARK_AS_ADVANCED(${_lib}_LIBRARY)

    IF(NOT ${_lib}_LIBRARY MATCHES "-NOTFOUND")
      LIST(APPEND BLAS_LIBRARIES ${${_lib}_LIBRARY})
      LIST(APPEND LAPACK_LIBRARIES ${${_lib}_LIBRARY})
    ENDIF()

  ENDFOREACH()

  #
  # Filter out spurious "FALSE" in the library lists:
  #
  IF(DEFINED BLAS_LIBRARIES)
    LIST(REMOVE_ITEM BLAS_LIBRARIES "FALSE")
  ENDIF()
  LIST(REMOVE_ITEM LAPACK_LIBRARIES "FALSE")

  MARK_AS_ADVANCED(
    BLAS_DIR
    LAPACK_DIR
    )

ELSE()
  SET(DEALII_LAPACK_FOUND FALSE)

  SET(LAPACK_DIR "" CACHE PATH
    "An optional hint to a LAPACK installation"
    )
  SET(BLAS_DIR "" CACHE PATH
    "An optional hint to a BLAS installation"
    )

  #
  # Clean up the library variables in case we couldn't find the libraries
  # to avoid spurious inclusions of "-NOTFOUND" or "FALSE":
  #
  SET(BLAS_LIBRARIES)
  SET(LAPACK_LIBRARIES)

ENDIF()
