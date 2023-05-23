## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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
#   TRILINOS_VERSION
#   TRILINOS_VERSION_MAJOR
#   TRILINOS_VERSION_MINOR
#   TRILINOS_VERSION_SUBMINOR
#   TRILINOS_WITH_MPI
#   TRILINOS_WITH_NO_64BITS_INDICES
#   TRILINOS_WITH_NO_32BITS_INDICES
#

set(TRILINOS_DIR "" CACHE PATH "An optional hint to a Trilinos installation")
set_if_empty(TRILINOS_DIR "$ENV{TRILINOS_DIR}")

#
# Include the trilinos package configuration:
#
find_package(TRILINOS_CONFIG
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


set(_target Trilinos::all_libs)
process_feature(TRILINOS
  TARGETS REQUIRED _target
  CLEAR
    TRILINOS_CONFIG_DIR EPETRA_CONFIG_H SACADO_CMATH_HPP ${_libraries}
    SACADO_CONFIG_H
    TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
  )


if(TRILINOS_FOUND)
  #
  # Extract version numbers:
  #
  set(TRILINOS_VERSION "${Trilinos_VERSION}")

  string(REGEX REPLACE
    "^([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MAJOR "${Trilinos_VERSION}")

  string(REGEX REPLACE
    "^[0-9]+\\.([0-9]+).*$" "\\1"
    TRILINOS_VERSION_MINOR "${Trilinos_VERSION}")

  # If there is no subminor number, TRILINOS_VERSION_SUBMINOR is set to an
  # empty string. If that is the case, set the subminor number to zero
  string(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.?(([0-9]+)?).*$" "\\1"
    TRILINOS_VERSION_SUBMINOR "${Trilinos_VERSION}")
  if("${TRILINOS_VERSION_SUBMINOR}" STREQUAL "")
    set(TRILINOS_VERSION_SUBMINOR "0")
  endif()

  #
  # Look for Epetra_config.h - we'll query it to determine MPI and 64bit
  # indices support:
  #
  deal_ii_find_file(EPETRA_CONFIG_H Epetra_config.h
    HINTS ${Trilinos_INCLUDE_DIRS}
    NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
    )

  if(EXISTS ${EPETRA_CONFIG_H})
    #
    # Determine whether Trilinos was configured with MPI and 64bit indices:
    #
    file(STRINGS "${EPETRA_CONFIG_H}" EPETRA_MPI_STRING
      REGEX "^[ \t]*#[ \t]*define[ \t]+HAVE_MPI")
    if("${EPETRA_MPI_STRING}" STREQUAL "")
      set(TRILINOS_WITH_MPI FALSE)
    else()
      set(TRILINOS_WITH_MPI TRUE)
    endif()
    file(STRINGS "${EPETRA_CONFIG_H}" EPETRA_32BIT_STRING
      REGEX "^[ \t]*#[ \t]*define[ \t]+EPETRA_NO_32BIT_GLOBAL_INDICES")
    if("${EPETRA_64BIT_STRING}" STREQUAL "")
      set(TRILINOS_WITH_NO_32BITS_INDICES TRUE)
    else()
      set(TRILINOS_WITH_NO_32BITS_INDICES FALSE)
    endif()
    file(STRINGS "${EPETRA_CONFIG_H}" EPETRA_64BIT_STRING
      REGEX "^[ \t]*#[ \t]*define[ \t]+EPETRA_NO_64BIT_GLOBAL_INDICES")
    if("${EPETRA_64BIT_STRING}" STREQUAL "")
      set(TRILINOS_WITH_NO_64BITS_INDICES TRUE)
    else()
      set(TRILINOS_WITH_NO_64BITS_INDICES FALSE)
    endif()
  endif()
endif()
