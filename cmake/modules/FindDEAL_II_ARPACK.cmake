## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# Try to find the ARPACK library
#
# This module exports
#
#   ARPACK_LIBRARIES
#   ARPACK_LINKER_FLAGS
#   ARPACK_WITH_PARPACK
#

set(ARPACK_DIR "" CACHE PATH "An optional hint to an ARPACK installation")
set_if_empty(ARPACK_DIR "$ENV{ARPACK_DIR}")

set(PARPACK_DIR "" CACHE PATH "An optional hint to a PARPACK installation")
set_if_empty(PARPACK_DIR "$ENV{PARPACK_DIR}")

deal_ii_find_library(ARPACK_LIBRARY
  NAMES arpack
  HINTS ${ARPACK_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

if(DEAL_II_WITH_MPI)
  get_filename_component(_path "${ARPACK_LIBRARY}" PATH)
  deal_ii_find_library(PARPACK_LIBRARY
    NAMES parpack
    HINTS ${_path} ${ARPACK_DIR} ${PARPACK_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
    )
else()
  set(PARPACK_LIBRARY "PARPACK_LIBRARY-NOTFOUND")
endif()

if(NOT DEAL_II_ARPACK_WITH_PARPACK)
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
  foreach(_libraries ${TRILINOS_LIBRARIES} ${PETSC_LIBRARIES})
    if("${_libraries}" MATCHES "scalapack")
      set(PARPACK_LIBRARY "PARPACK_LIBRARY-NOTFOUND")
    endif()
  endforeach()
endif()


if(NOT PARPACK_LIBRARY MATCHES "-NOTFOUND")
  set(ARPACK_WITH_PARPACK TRUE)
else()
  set(ARPACK_WITH_PARPACK FALSE)
endif()

process_feature(ARPACK
  LIBRARIES
    OPTIONAL PARPACK_LIBRARY
    REQUIRED ARPACK_LIBRARY LAPACK_LIBRARIES
    OPTIONAL MPI_C_LIBRARIES
  LINKER_FLAGS OPTIONAL LAPACK_LINKER_FLAGS
  CLEAR ARPACK_LIBRARY PARPACK_LIBRARY
  )
