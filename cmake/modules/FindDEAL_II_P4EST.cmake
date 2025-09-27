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
# Try to find the P4EST library
#
# This module exports:
#   P4EST_LIBRARIES
#   P4EST_INCLUDE_DIRS
#   P4EST_WITH_MPI
#   P4EST_WITH_ZLIB
#   P4EST_WITH_VTK_BINARY
#   P4EST_VERSION
#   P4EST_VERSION_MAJOR
#   P4EST_VERSION_MINOR
#   P4EST_VERSION_SUBMINOR
#   P4EST_VERSION_PATCH
#

set(P4EST_DIR "" CACHE PATH
  "An optional hint to a p4est installation/directory"
  )
set_if_empty(P4EST_DIR "$ENV{P4EST_DIR}")
set_if_empty(SC_DIR "$ENV{SC_DIR}")

#
# Search for the sc library, usually bundled with p4est. If no SC_DIR was
# given, take what we chose for p4est.
#

deal_ii_find_path(SC_INCLUDE_DIR sc.h
  HINTS
    ${SC_DIR}/FAST
    ${SC_DIR}/DEBUG
    ${SC_DIR}
    ${P4EST_DIR}/FAST
    ${P4EST_DIR}/DEBUG
    ${P4EST_DIR}
  PATH_SUFFIXES
    sc include/p4est include src sc/src
  )

deal_ii_find_library(P4EST_LIBRARY_OPTIMIZED
  NAMES p4est
  HINTS ${P4EST_DIR}/FAST ${P4EST_DIR}/DEBUG ${P4EST_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib src
  )

deal_ii_find_library(SC_LIBRARY_OPTIMIZED
  NAMES sc
  HINTS
    ${SC_DIR}/FAST
    ${SC_DIR}/DEBUG
    ${SC_DIR}
    ${P4EST_DIR}/FAST
    ${P4EST_DIR}/DEBUG
    ${P4EST_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib src sc/src
  )

#
# Support debug variants as well:
#

deal_ii_find_library(P4EST_LIBRARY_DEBUG
  NAMES p4est
  HINTS ${P4EST_DIR}/DEBUG
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib src
  )

deal_ii_find_library(SC_LIBRARY_DEBUG
  NAMES sc
  HINTS ${SC_DIR}/DEBUG ${P4EST_DIR}/DEBUG
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib src sc/src
  )

if( ( "${P4EST_LIBRARY_OPTIMIZED}" STREQUAL "${P4EST_LIBRARY_DEBUG}"
      AND "${SC_LIBRARY_OPTIMIZED}" STREQUAL "${SC_LIBRARY_DEBUG}" )
    OR P4EST_LIBRARY_DEBUG MATCHES "-NOTFOUND"
    OR SC_LIBRARY_DEBUG MATCHES "-NOTFOUND" )
  set(_libraries
    LIBRARIES REQUIRED P4EST_LIBRARY_OPTIMIZED SC_LIBRARY_OPTIMIZED
    )
else()
  set(_libraries
    LIBRARIES_RELEASE REQUIRED P4EST_LIBRARY_OPTIMIZED SC_LIBRARY_OPTIMIZED
    LIBRARIES_DEBUG REQUIRED P4EST_LIBRARY_DEBUG SC_LIBRARY_DEBUG
    )
endif()


deal_ii_find_path(P4EST_INCLUDE_DIR p4est_config.h
  HINTS ${P4EST_DIR}/FAST ${P4EST_DIR}/DEBUG ${P4EST_DIR}
  PATH_SUFFIXES p4est include/p4est include src
  )

if(EXISTS ${P4EST_INCLUDE_DIR}/p4est_config.h)
  #
  # Determine mpi support of p4est:
  #
  file(STRINGS "${P4EST_INCLUDE_DIR}/p4est_config.h" P4EST_MPI_STRING
    REGEX "#define.*P4EST_MPI 1")
  if("${P4EST_MPI_STRING}" STREQUAL "")
    file(STRINGS "${P4EST_INCLUDE_DIR}/p4est_config.h" P4EST_MPI_STRING
      REGEX "#define.*P4EST_ENABLE_MPI")
    if("${P4EST_MPI_STRING}" STREQUAL "")
      set(P4EST_WITH_MPI FALSE)
    else()
      set(P4EST_WITH_MPI TRUE)
    endif()
  else()
    set(P4EST_WITH_MPI TRUE)
  endif()

  #
  # Is p4est built against zlib?
  #
  file(STRINGS "${P4EST_INCLUDE_DIR}/p4est_config.h" P4EST_ZLIB_STRING
    REGEX "^#define.*P4EST_HAVE_ZLIB")
  if("${P4EST_ZLIB_STRING}" STREQUAL "")
    set(P4EST_WITH_ZLIB FALSE)
  else()
    set(P4EST_WITH_ZLIB TRUE)
  endif()

  #
  # Is binary vtk output enabled?
  #
  file(STRINGS "${P4EST_INCLUDE_DIR}/p4est_config.h" P4EST_VTK_BINARY_STRING
    REGEX "#define.*P4EST_ENABLE_VTK_BINARY 1")
  if("${P4EST_VTK_BINARY_STRING}" STREQUAL "")
    set(P4EST_WITH_VTK_BINARY FALSE)
  else()
    set(P4EST_WITH_VTK_BINARY TRUE)
  endif()

  #
  # Extract version numbers:
  #
  file(STRINGS "${P4EST_INCLUDE_DIR}/p4est_config.h" P4EST_VERSION
    REGEX "^[ \t]*#[ \t]*define[ \t]+P4EST_VERSION \"")
  string(REGEX REPLACE "^.*P4EST_VERSION.*\"([0-9]+.*)\".*" "\\1"
    P4EST_VERSION "${P4EST_VERSION}"
    )
  string(REGEX REPLACE
    "^([0-9]+).*$" "\\1"
    P4EST_VERSION_MAJOR "${P4EST_VERSION}")
  string(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*$" "\\1"
    P4EST_VERSION_MINOR "${P4EST_VERSION}")
  string(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
    P4EST_VERSION_SUBMINOR "${P4EST_VERSION}")
  string(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
    P4EST_VERSION_PATCH "${P4EST_VERSION}")

  #
  # We cannot rely on the fact that SUBMINOR or PATCH are defined.
  # Nevertheless, we need a full version number for our preprocessor macros
  # to work. If the p4est version number is only of the form x.y instead of
  # a.b.c.d, then the last two REGEX_REPLACE calls above will have failed
  # because the regular expression didn't match the version string,
  # and P4EST_VERSION_SUBMINOR and P4EST_VERSION_PATCH will either be
  # empty or be the full version string. In those cases, set those numbers
  # to 0 if necessary.
  #
  if("${P4EST_VERSION_SUBMINOR}" MATCHES "^(|${P4EST_VERSION})$")
    set(P4EST_VERSION_SUBMINOR "0")
  endif()

  if("${P4EST_VERSION_PATCH}" MATCHES "^(|${P4EST_VERSION})$")
    set(P4EST_VERSION_PATCH "0")
  endif()
endif()

process_feature(P4EST
  ${_libraries}
  LIBRARIES OPTIONAL LAPACK_LIBRARIES MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED P4EST_INCLUDE_DIR SC_INCLUDE_DIR
  CLEAR
    SC_INCLUDE_DIR P4EST_LIBRARY_OPTIMIZED SC_LIBRARY_OPTIMIZED
    P4EST_LIBRARY_DEBUG SC_LIBRARY_DEBUG P4EST_INCLUDE_DIR
  )
