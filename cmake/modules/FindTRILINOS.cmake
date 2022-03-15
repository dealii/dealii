## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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
# Try to find the Trilinos library
#
# This module exports:
#
#   TRILINOS_DIR
#   TRILINOS_INCLUDE_DIRS
#   TRILINOS_LIBRARIES
#   TRILINOS_LINKER_FLAGS
#   TRILINOS_VERSION
#   TRILINOS_VERSION_MAJOR
#   TRILINOS_VERSION_MINOR
#   TRILINOS_VERSION_SUBMINOR
#   TRILINOS_WITH_MPI
#

SET(TRILINOS_DIR "" CACHE PATH "An optional hint to a Trilinos installation")
SET_IF_EMPTY(TRILINOS_DIR "$ENV{TRILINOS_DIR}")

#
# Include the trilinos package configuration:
#
FIND_PACKAGE(TRILINOS_CONFIG
  CONFIG QUIET
  NAMES Trilinos TRILINOS
  HINTS
    ${TRILINOS_DIR}/lib/cmake/Trilinos
    ${TRILINOS_DIR}
  PATH_SUFFIXES
    lib64/cmake/Trilinos
    lib/cmake/Trilinos
    lib${LIB_SUFFIX}/cmake/Trilinos
  NO_SYSTEM_ENVIRONMENT_PATH
  )


IF(DEFINED Trilinos_VERSION)
  #
  # Extract version numbers:
  #
  SET(TRILINOS_VERSION "${Trilinos_VERSION}")

  STRING(REGEX REPLACE
    "^([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MAJOR "${Trilinos_VERSION}")

  STRING(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MINOR "${Trilinos_VERSION}")

  # If there is no subminor number, 
  # TRILINOS_VERSION_SUBMINOR is set to an empty string. 
  # If that is the case, set the subminor number to zero
  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.?(([0-9]+)?).*$" "\\1"
    TRILINOS_VERSION_SUBMINOR "${Trilinos_VERSION}")
  IF("${TRILINOS_VERSION_SUBMINOR}" STREQUAL "")
    SET(TRILINOS_VERSION_SUBMINOR "0")
  ENDIF()  
ENDIF()

#
# Look for Epetra_config.h - we'll query it to determine MPI and 64bit
# indices support:
#
DEAL_II_FIND_FILE(EPETRA_CONFIG_H Epetra_config.h
  HINTS ${Trilinos_INCLUDE_DIRS}
  NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
  )

IF(EXISTS ${EPETRA_CONFIG_H})
  #
  # Determine whether Trilinos was configured with MPI and 64bit indices:
  #
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_MPI_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+HAVE_MPI")
  IF("${EPETRA_MPI_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_MPI FALSE)
  ELSE()
    SET(TRILINOS_WITH_MPI TRUE)
  ENDIF()
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_32BIT_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+EPETRA_NO_32BIT_GLOBAL_INDICES")
  IF("${EPETRA_64BIT_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_NO_32BITS_INDICES TRUE)
  ELSE()
    SET(TRILINOS_WITH_NO_32BITS_INDICES FALSE)
  ENDIF()
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_64BIT_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+EPETRA_NO_64BIT_GLOBAL_INDICES")
  IF("${EPETRA_64BIT_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_NO_64BITS_INDICES TRUE)
  ELSE()
    SET(TRILINOS_WITH_NO_64BITS_INDICES FALSE)
  ENDIF()
ENDIF()


#
# *Boy* Sanitize variables that are exported by TrilinosConfig.cmake...
#
# Especially deduplicate stuff...
#
REMOVE_DUPLICATES(Trilinos_LIBRARIES REVERSE)
REMOVE_DUPLICATES(Trilinos_TPL_LIBRARIES REVERSE)

REMOVE_DUPLICATES(Trilinos_INCLUDE_DIRS)
STRING(REGEX REPLACE
  "(lib64|lib)\\/cmake\\/Trilinos\\/\\.\\.\\/\\.\\.\\/\\.\\.\\/" ""
  Trilinos_INCLUDE_DIRS "${Trilinos_INCLUDE_DIRS}"
  )

REMOVE_DUPLICATES(Trilinos_TPL_INCLUDE_DIRS)

#
# We'd like to have the full library names but the Trilinos package only
# exports a list with short names...
# So we check again for every lib and store the full path:
#
SET(_libraries "")
FOREACH(_library ${Trilinos_LIBRARIES})
  LIST(APPEND _libraries TRILINOS_LIBRARY_${_library})
  DEAL_II_FIND_LIBRARY(TRILINOS_LIBRARY_${_library}
    NAMES ${_library}
    HINTS ${Trilinos_LIBRARY_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    NO_CMAKE_FIND_ROOT_PATH
    )
ENDFOREACH()


DEAL_II_PACKAGE_HANDLE(TRILINOS
  LIBRARIES
    REQUIRED ${_libraries}
    OPTIONAL Trilinos_TPL_LIBRARIES MPI_CXX_LIBRARIES
  INCLUDE_DIRS
    REQUIRED Trilinos_INCLUDE_DIRS
    OPTIONAL Trilinos_TPL_INCLUDE_DIRS
  USER_INCLUDE_DIRS
    REQUIRED Trilinos_INCLUDE_DIRS
    OPTIONAL Trilinos_TPL_INCLUDE_DIRS
  LINKER_FLAGS
    OPTIONAL Trilinos_EXTRA_LD_FLAGS
  CLEAR
    TRILINOS_CONFIG_DIR EPETRA_CONFIG_H SACADO_CMATH_HPP ${_libraries}
    SACADO_CONFIG_H
    TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
  )
