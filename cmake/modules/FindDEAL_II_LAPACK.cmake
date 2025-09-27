## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
option(LAPACK_WITH_64BIT_BLAS_INDICES
  "BLAS has 64 bit integers."
  OFF
  )
mark_as_advanced(LAPACK_WITH_64BIT_BLAS_INDICES)

set(LAPACK_DIR "" CACHE PATH "An optional hint to a LAPACK installation")
set(BLAS_DIR "" CACHE PATH "An optional hint to a BLAS installation")
set_if_empty(BLAS_DIR "$ENV{BLAS_DIR}")
set_if_empty(LAPACK_DIR "$ENV{LAPACK_DIR}")

set(_cmake_prefix_path_backup "${CMAKE_PREFIX_PATH}")

set(CMAKE_PREFIX_PATH ${BLAS_DIR} ${LAPACK_DIR} ${_cmake_prefix_path_backup})
find_package(BLAS)

set(CMAKE_PREFIX_PATH ${LAPACK_DIR} ${_cmake_prefix_path_backup})
find_package(LAPACK)

set(CMAKE_PREFIX_PATH ${_cmake_prefix_path_backup})

#
# Filter out spurious "FALSE" in the library lists:
#
if(DEFINED BLAS_LIBRARIES)
  list(REMOVE_ITEM BLAS_LIBRARIES "FALSE")
endif()
if(DEFINED LAPACK_LIBRARIES)
  list(REMOVE_ITEM LAPACK_LIBRARIES "FALSE")
endif()

#
# Work around a bug in CMake 3.11 by simply filtering out
# "PkgConf::PKGC_BLAS". See bug
#   https://gitlab.kitware.com/cmake/cmake/issues/17934
#
if(DEFINED BLAS_LIBRARIES)
  list(REMOVE_ITEM BLAS_LIBRARIES "PkgConfig::PKGC_BLAS")
endif()

#
# Well, in case of static archives we have to manually pick up the
# complete link interface. *sigh*
#
# If CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES is not available, do it
# unconditionally for the most common case (gfortran).
#
if(NOT BUILD_SHARED_LIBS)
  set(_fortran_libs ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  #
  # Since CMake 3.9 the gcc runtime libraries libgcc.a and libgcc_s.so.1
  # have been added to the CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES variable.
  # We thus have to remove the shared low-level runtime library
  # libgcc_s.so.1 from the link interface; otherwise completely static
  # linkage is broken.
  #
  list(REMOVE_ITEM _fortran_libs gcc_s)
  set_if_empty(_fortran_libs gfortran quadmath m)

  foreach(_lib ${_fortran_libs})
    find_system_library(${_lib}_LIBRARY NAMES ${_lib})
    list(APPEND _additional_libraries ${_lib}_LIBRARY)
  endforeach()
endif()


process_feature(LAPACK
  LIBRARIES
    REQUIRED LAPACK_LIBRARIES
    OPTIONAL BLAS_LIBRARIES ${_additional_libraries}
  LINKER_FLAGS OPTIONAL LAPACK_LINKER_FLAGS BLAS_LINKER_FLAGS
  INCLUDE_DIRS OPTIONAL LAPACK_INCLUDE_DIRS
  CLEAR
    atlas_LIBRARY atlcblas_LIBRARY atllapack_LIBRARY blas_LIBRARY
    eigen_blas_LIBRARY f77blas_LIBRARY gslcblas_LIBRARY lapack_LIBRARY
    m_LIBRARY ptf77blas_LIBRARY ptlapack_LIBRARY refblas_LIBRARY
    reflapack_LIBRARY BLAS_LIBRARIES ${_additional_libraries}
    LAPACK_SYMBOL_CHECK # clean up check in configure_1_lapack.cmake
  )
