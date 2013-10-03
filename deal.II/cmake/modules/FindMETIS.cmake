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

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(METIS_DIR "$ENV{METIS_DIR}")

#
# Metis is usually pretty self contained. So no external dependencies
# so far. But there could be dependencies on pcre and mpi...
#
# Link in MPI unconditionally (if found).
#

FIND_PATH(METIS_INCLUDE_DIR metis.h
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    metis include/metis include
  )

FIND_LIBRARY(METIS_LIBRARY
  NAMES metis
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
  )

FIND_LIBRARY(PARMETIS_LIBRARY
  NAMES parmetis
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(METIS DEFAULT_MSG
  METIS_LIBRARY
  METIS_INCLUDE_DIR
  )

MARK_AS_ADVANCED(
  METIS_LIBRARY
  PARMETIS_LIBRARY
  METIS_INCLUDE_DIR
  )

IF(METIS_FOUND)

  #
  # Sanity check: Only include parmetis library if it is in the same
  # directory as the metis library...
  #
  GET_FILENAME_COMPONENT(_path1 "${PARMETIS_LIBRARY}" PATH)
  GET_FILENAME_COMPONENT(_path2 "${ETIS_LIBRARY}" PATH)
  IF( NOT PARMETIS_LIBRARY MATCHES "-NOTFOUND"
      AND "${_path1}" STREQUAL "${_path2}" )
    SET(METIS_LIBRARIES ${PARMETIS_LIBRARY})
  ENDIF()

  LIST(APPEND METIS_LIBRARIES
    ${METIS_LIBRARY}
    ${MPI_C_LIBRARIES} # for good measure
    )
  SET(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})

  #
  # Extract the version number out of metis.h
  #
  FILE(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_major_string
    REGEX "METIS_VER_MAJOR"
    )
  STRING(REGEX REPLACE "^.*METIS_VER_MAJOR.* ([0-9]+).*" "\\1"
    METIS_VERSION_MAJOR "${_metis_major_string}"
    )
  FILE(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_minor_string
    REGEX "METIS_VER_MINOR"
    )
  STRING(REGEX REPLACE "^.*METIS_VER_MINOR.* ([0-9]+).*" "\\1"
    METIS_VERSION_MINOR "${_metis_minor_string}"
    )
  FILE(STRINGS "${METIS_INCLUDE_DIR}/metis.h" _metis_subminor_string
    REGEX "METIS_VER_SUBMINOR"
    )
  STRING(REGEX REPLACE "^.*METIS_VER_SUBMINOR.* ([0-9]+).*" "\\1"
    METIS_VERSION_SUBMINOR "${_metis_subminor_string}"
    )
  SET(METIS_VERSION
    "${METIS_VERSION_MAJOR}.${METIS_VERSION_MINOR}.${METIS_VERSION_SUBMINOR}"
    )

  MARK_AS_ADVANCED(METIS_DIR)
ELSE()
  SET(METIS_DIR "" CACHE PATH
    "An optional hint to a metis directory"
    )
ENDIF()

