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

if(Trilinos_FOUND)
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

if(TRILINOS_VERSION VERSION_GREATER_EQUAL 14)
  set(_targets Trilinos::all_libs)
  if(TARGET Kokkos::kokkos)
    list(APPEND _targets Kokkos::kokkos)
  endif()

  process_feature(TRILINOS
    TARGETS REQUIRED _targets
    CLEAR
      EPETRA_CONFIG_H TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
    )
else()
  #
  # *Boy* Sanitize variables that are exported by TrilinosConfig.cmake...
  #
  # Especially deduplicate stuff...
  #
  remove_duplicates(Trilinos_LIBRARIES REVERSE)
  remove_duplicates(Trilinos_TPL_LIBRARIES REVERSE)

  remove_duplicates(Trilinos_INCLUDE_DIRS)
  string(REGEX REPLACE
    "(lib64|lib)\\/cmake\\/Trilinos\\/\\.\\.\\/\\.\\.\\/\\.\\.\\/" ""
    Trilinos_INCLUDE_DIRS "${Trilinos_INCLUDE_DIRS}"
    )

  remove_duplicates(Trilinos_TPL_INCLUDE_DIRS)

  if(TARGET Kokkos::kokkos)
    get_property(KOKKOS_COMPILE_FLAGS_FULL TARGET Kokkos::kokkos PROPERTY INTERFACE_COMPILE_OPTIONS)
    string(REGEX REPLACE "\\$<\\$<COMPILE_LANGUAGE:CXX>:([^>]*)>" "\\1" KOKKOS_COMPILE_FLAGS "${KOKKOS_COMPILE_FLAGS_FULL}")
    string(REPLACE ";" " " KOKKOS_COMPILE_FLAGS "${KOKKOS_COMPILE_FLAGS}")

    #
    # Extract missing openmp linker options from the
    # "INTERFACE_LINK_LIBRARIES" property of the Kokkos::kokkos target:
    #
    set(_kokkos_openmp_flags)
    get_property(_libraries TARGET Kokkos::kokkos PROPERTY INTERFACE_LINK_LIBRARIES)
    foreach(_entry ${_libraries})
      if(${_entry} MATCHES "^-f?openmp")
        add_flags(_kokkos_openmp_flags "${_entry}")
      endif()
    endforeach()

    # Some Kokkos versions included in Trilinos before 13.2.0 add "-x cuda" when 
    # using clang++ as compiler even if Kokkos has not been configured with Cuda
    # support. Simply strip that flag from what we are using in that case.
    if(NOT Kokkos_ENABLE_CUDA)
      string(REPLACE "-x cuda" "" KOKKOS_COMPILE_FLAGS "${KOKKOS_COMPILE_FLAGS}")
    endif()
  endif()

  #
  # We'd like to have the full library names but the Trilinos package only
  # exports a list with short names...
  # So we check again for every lib and store the full path:
  #
  set(_libraries "")
  foreach(_library ${Trilinos_LIBRARIES})
    list(APPEND _libraries TRILINOS_LIBRARY_${_library})
    deal_ii_find_library(TRILINOS_LIBRARY_${_library}
      NAMES ${_library}
      HINTS ${Trilinos_LIBRARY_DIRS}
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH
      )
  endforeach()

  process_feature(TRILINOS
    LIBRARIES
      REQUIRED ${_libraries}
      OPTIONAL Trilinos_TPL_LIBRARIES MPI_CXX_LIBRARIES
    INCLUDE_DIRS
      REQUIRED Trilinos_INCLUDE_DIRS
      OPTIONAL Trilinos_TPL_INCLUDE_DIRS
    CXX_FLAGS
      OPTIONAL KOKKOS_COMPILE_FLAGS
    LINKER_FLAGS
    OPTIONAL Trilinos_EXTRA_LD_FLAGS _kokkos_openmp_flags
    CLEAR
      EPETRA_CONFIG_H ${_libraries} TRILINOS_CXX_SUPPORTS_SACADO_COMPLEX_RAD
    )
endif()
