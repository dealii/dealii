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
# Try to find the (serial) METIS library
#
# This module exports
#
#   METIS_LIBRARIES
#   METIS_INCLUDE_DIRS
#   METIS_VERSION
#   METIS_VERSION_MAJOR
#   METIS_VERSION_MINOR
#   METIS_VERSION_SUBMINOR
#

set(METIS_DIR "" CACHE PATH "An optional hint to a metis directory")
set_if_empty(METIS_DIR "$ENV{METIS_DIR}")

#
# Metis is usually pretty self contained. So no external dependencies
# so far. But there could be dependencies on pcre and mpi...
#
# Link in MPI unconditionally (if found).
#

deal_ii_find_library(METIS_LIBRARY
  NAMES metis
  HINTS ${METIS_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
  )


deal_ii_find_path(METIS_INCLUDE_DIR metis.h
  HINTS ${METIS_DIR}
  PATH_SUFFIXES metis include/metis include
  )

if(EXISTS ${METIS_INCLUDE_DIR}/metis.h)
  #
  # Extract the version number out of metis.h
  #
  file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_major_string
    REGEX "METIS_VER_MAJOR"
    )
  string(REGEX REPLACE "^.*METIS_VER_MAJOR.* ([0-9]+).*" "\\1"
    METIS_VERSION_MAJOR "${_metis_major_string}"
    )
  file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_minor_string
    REGEX "METIS_VER_MINOR"
    )
  string(REGEX REPLACE "^.*METIS_VER_MINOR.* ([0-9]+).*" "\\1"
    METIS_VERSION_MINOR "${_metis_minor_string}"
    )
  file(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_subminor_string
    REGEX "METIS_VER_SUBMINOR"
    )
  string(REGEX REPLACE "^.*METIS_VER_SUBMINOR.* ([0-9]+).*" "\\1"
    METIS_VERSION_SUBMINOR "${_metis_subminor_string}"
    )
  set(METIS_VERSION
    "${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_SUBMINOR}"
    )
  if("${METIS_VERSION}" STREQUAL "..")
    set(METIS_VERSION)
  endif()
endif()

process_feature(METIS
  LIBRARIES
    REQUIRED METIS_LIBRARY
    OPTIONAL MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED METIS_INCLUDE_DIR
  CLEAR METIS_LIBRARY METIS_INCLUDE_DIR
  )
