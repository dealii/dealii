## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2014 by the deal.II authors
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
SET(LAPACK_DIR "" CACHE PATH "An optional hint to a LAPACK installation")
SET(BLAS_DIR "" CACHE PATH "An optional hint to a BLAS installation")
SET_IF_EMPTY(BLAS_DIR "$ENV{BLAS_DIR}")
SET_IF_EMPTY(LAPACK_DIR "$ENV{LAPACK_DIR}")

SET(_cmake_prefix_path_backup "${CMAKE_PREFIX_PATH}")
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

SET(CMAKE_PREFIX_PATH ${BLAS_DIR} ${LAPACK_DIR} ${_cmake_prefix_path_backup})
FIND_PACKAGE(BLAS)

SET(CMAKE_PREFIX_PATH ${LAPACK_DIR} ${_cmake_prefix_path_backup})
FIND_PACKAGE(LAPACK)

SET(CMAKE_PREFIX_PATH ${_cmake_prefix_path_backup})
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

#
# Filter out spurious "FALSE" in the library lists:
#
IF(DEFINED BLAS_LIBRARIES)
  LIST(REMOVE_ITEM BLAS_LIBRARIES "FALSE")
ENDIF()
IF(DEFINED LAPACK_LIBRARIES)
  LIST(REMOVE_ITEM LAPACK_LIBRARIES "FALSE")
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
  LIST(APPEND _additional_libraries ${_lib}_LIBRARY)
ENDFOREACH()


SET(_lapack_libraries ${LAPACK_LIBRARIES})
SET(_lapack_linker_flags ${LAPACK_LINKER_FLAGS})
DEAL_II_PACKAGE_HANDLE(LAPACK
  LIBRARIES
    REQUIRED _lapack_libraries
    OPTIONAL BLAS_LIBRARIES ${_additional_libraries}
  LINKER_FLAGS OPTIONAL _lapack_linker_flags BLAS_LINKER_FLAGS
  CLEAR
    atlas_LIBRARY atlcblas_LIBRARY atllapack_LIBRARY blas_LIBRARY
    eigen_blas_LIBRARY f77blas_LIBRARY gslcblas_LIBRARY lapack_LIBRARY
    m_LIBRARY ptf77blas_LIBRARY ptlapack_LIBRARY refblas_LIBRARY
    reflapack_LIBRARY BLAS_LIBRARIES ${_additional_libraries}
    LAPACK_SYMBOL_CHECK # Cleanup check in configure_1_lapack.cmake
  )
