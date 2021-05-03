## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2019 by the deal.II authors
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
OPTION(LAPACK_WITH_64BIT_BLAS_INDICES
  "BLAS has 64 bit integers."
  OFF
  )
MARK_AS_ADVANCED(LAPACK_WITH_64BIT_BLAS_INDICES)

SET(LAPACK_DIR "" CACHE PATH "An optional hint to a LAPACK installation")
SET(BLAS_DIR "" CACHE PATH "An optional hint to a BLAS installation")
SET_IF_EMPTY(BLAS_DIR "$ENV{BLAS_DIR}")
SET_IF_EMPTY(LAPACK_DIR "$ENV{LAPACK_DIR}")

SET(_cmake_prefix_path_backup "${CMAKE_PREFIX_PATH}")

# temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
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
# Work around a bug in CMake 3.11 by simply filtering out
# "PkgConf::PKGC_BLAS". See bug
#   https://gitlab.kitware.com/cmake/cmake/issues/17934
#
IF(DEFINED BLAS_LIBRARIES)
  LIST(REMOVE_ITEM BLAS_LIBRARIES "PkgConfig::PKGC_BLAS")
ENDIF()

#
# Well, in case of static archives we have to manually pick up the
# complete link interface. *sigh*
#
# If CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is not available, do it
# unconditionally for the most common case (gfortran).
#
IF(NOT BUILD_SHARED_LIBS)
  SET(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  #
  # Since CMake 3.9 the gcc runtime libraries libgcc.a and libgcc_s.so.1
  # have been added to the CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES variable.
  # We thus have to remove the shared low-level runtime library
  # libgcc_s.so.1 from the link interface; otherwise completely static
  # linkage is broken.
  #
  LIST(REMOVE_ITEM _fortran_libs gcc_s)
  SET_IF_EMPTY(_fortran_libs gfortran quadmath m)

  FOREACH(_lib ${_fortran_libs})
    FIND_SYSTEM_LIBRARY(${_lib}_LIBRARY NAMES ${_lib})
    LIST(APPEND _additional_libraries ${_lib}_LIBRARY)
  ENDFOREACH()
ENDIF()


SET(_lapack_include_dirs ${LAPACK_INCLUDE_DIRS})
SET(_lapack_libraries ${LAPACK_LIBRARIES})
SET(_lapack_linker_flags ${LAPACK_LINKER_FLAGS})
DEAL_II_PACKAGE_HANDLE(LAPACK
  LIBRARIES
    REQUIRED _lapack_libraries
    OPTIONAL BLAS_LIBRARIES ${_additional_libraries}
  LINKER_FLAGS OPTIONAL _lapack_linker_flags BLAS_LINKER_FLAGS
  INCLUDE_DIRS
    OPTIONAL _lapack_include_dirs
  USER_INCLUDE_DIRS
    OPTIONAL _lapack_include_dirs
  CLEAR
    atlas_LIBRARY atlcblas_LIBRARY atllapack_LIBRARY blas_LIBRARY
    eigen_blas_LIBRARY f77blas_LIBRARY gslcblas_LIBRARY lapack_LIBRARY
    m_LIBRARY ptf77blas_LIBRARY ptlapack_LIBRARY refblas_LIBRARY
    reflapack_LIBRARY BLAS_LIBRARIES ${_additional_libraries}
    LAPACK_SYMBOL_CHECK # clean up check in configure_1_lapack.cmake
  )
