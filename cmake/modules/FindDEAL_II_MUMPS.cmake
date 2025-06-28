## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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
# Try to find the MUMPS library
#
# This module exports
#
#   MUMPS_INCLUDE_DIRS
#   MUMPS_LIBRARIES
#   MUMPS_LINKER_FLAGS
#   MUMPS_VERSION
#   MUMPS_VERSION_MAJOR
#   MUMPS_VERSION_MINOR
#   MUMPS_VERSION_SUBMINOR
#

set(MUMPS_DIR "" CACHE PATH "An optional hint to a mumps directory")
set_if_empty(MUMPS_DIR "$ENV{MUMPS_DIR}")

#
# Search for mumps:
#

deal_ii_find_path(MUMPS_INCLUDE_DIR dmumps_c.h
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES mumps include/mumps include
  )

deal_ii_find_library(DMUMPS_LIBRARY
  NAMES dmumps
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_library(MUMPS_COMMON_LIBRARY
  NAMES mumps_common
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# If we can find libpord.so (or similar), link it in as well:
#
deal_ii_find_library(PORD_LIBRARY
  NAMES pord
  HINTS ${MUMPS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

if(EXISTS ${MUMPS_INCLUDE_DIR}/dmumps_c.h)
  file(STRINGS "${MUMPS_INCLUDE_DIR}/dmumps_c.h" MUMPS_VERSION_STRING
    REGEX "#define.*MUMPS_VERSION")
  string(REGEX REPLACE "^.*MUMPS_VERSION.*\"(.+)\".*" "\\1"
    MUMPS_VERSION "${MUMPS_VERSION_STRING}"
    )
  string(REGEX REPLACE
    "([0-9]+)\\..*" "\\1" MUMPS_VERSION_MAJOR "${MUMPS_VERSION}"
    )
  string(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*" "\\1" MUMPS_VERSION_MINOR "${MUMPS_VERSION}"
    )
  string(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" MUMPS_VERSION_SUBMINOR "${MUMPS_VERSION}"
    )
endif()

process_feature(MUMPS
  LIBRARIES
    REQUIRED DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY SCALAPACK_LIBRARIES
    OPTIONAL BLACS_LIBRARIES
    OPTIONAL PORD_LIBRARY
    OPTIONAL METIS_LIBRARIES MPI_Fortran_LIBRARIES
  INCLUDE_DIRS
    REQUIRED MUMPS_INCLUDE_DIR
  LINKER_FLAGS
    OPTIONAL LAPACK_LINKER_FLAGS
  CLEAR
    DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY PORD_LIBRARY
    MUMPS_INCLUDE_DIR
  )
