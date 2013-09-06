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
# Try to find the Trilinos library
#
# This module exports:
#
#   TRILINOS_DIR
#   TRILINOS_INCLUDE_DIRS
#   TRILINOS_LIBRARIES
#   TRILINOS_VERSION
#   TRILINOS_VERSION_MAJOR
#   TRILINOS_VERSION_MINOR
#   TRILINOS_VERSION_SUBMINOR
#   TRILINOS_WITH_MPI
#   TRILINOS_SUPPORTS_CPP11
#   TRILINOS_HAS_C99_TR1_WORKAROUND
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(TRILINOS_DIR "$ENV{TRILINOS_DIR}")

#
# Do not include TrilinosConfig.cmake directly, it is just too big o_O
#
# Just search for the file:
#
FIND_FILE(TRILINOS_CONFIG
  NAMES TrilinosConfig.cmake trilinos-config.cmake
  HINTS
    ${TRILINOS_DIR}
  PATH_SUFFIXES
    lib64/cmake/Trilinos
    lib/cmake/Trilinos
    lib${LIB_SUFFIX}/cmake/Trilinos
    include/trilinos
    include/Trilinos
  NO_SYSTEM_ENVIRONMENT_PATH
  )

IF(NOT "${TRILINOS_CONFIG}" STREQUAL "${TRILINOS_CONFIG_SAVED}")
  SET(_new_trilinos_config TRUE)
ENDIF()
SET(TRILINOS_CONFIG_SAVED "${TRILINOS_CONFIG}" CACHE INTERNAL "" FORCE)


IF(NOT TRILINOS_CONFIG MATCHES "-NOTFOUND")

  SET(_filtered_trilinos_config
    "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/TrilinosConfig.cmake"
    )

  IF(_new_trilinos_config)
    GET_FILENAME_COMPONENT(_trilinos_path "${TRILINOS_CONFIG}" PATH)
    FILE(WRITE ${_filtered_trilinos_config} "SET(_cmake_current_list_dir ${_trilinos_path})\n")

    #
    # Only pick up every line that starts with "^SET("...
    #
    FILE(STRINGS "${TRILINOS_CONFIG}" _trilinos_config_filtered REGEX "^SET")

    FOREACH(_line ${_trilinos_config_filtered})
      STRING(REPLACE "CMAKE_CURRENT_LIST_DIR" "_cmake_current_list_dir"
        _line "${_line}"
        )
      FILE(APPEND ${_filtered_trilinos_config} "${_line}\n")
    ENDFOREACH()
  ENDIF()

  #
  # ... and include only that:
  #
  INCLUDE(${_filtered_trilinos_config})

  SET(TRILINOS_CONFIG_FOUND TRUE)
ENDIF()


#
# Look for the one include file that we'll query for further information:
#
IF(_new_trilinos_config)
  UNSET(EPETRA_CONFIG_H CACHE)
ENDIF()
FIND_FILE(EPETRA_CONFIG_H Epetra_config.h
  HINTS ${Trilinos_INCLUDE_DIRS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  NO_CMAKE_FIND_ROOT_PATH
  )
IF(EPETRA_CONFIG_H MATCHES "-NOTFOUND")
  SET(TRILINOS_CONFIG_FOUND FALSE)
ELSE()
  SET(TRILINOS_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS})
  #
  # *Boy* Sanitize the include paths given by TrilinosConfig.cmake...
  #
  STRING(REGEX REPLACE
    "(lib64|lib)\\/cmake\\/Trilinos\\/\\.\\.\\/\\.\\.\\/\\.\\.\\/" ""
    TRILINOS_INCLUDE_DIRS "${TRILINOS_INCLUDE_DIRS}"
    )
ENDIF()


#
# We'd like to have the full library names but the Trilinos package only
# exports a list with short names...
# So we check again for every lib and store the full path:
#
FOREACH(_library ${Trilinos_LIBRARIES})
  IF(_new_trilinos_config)
    UNSET(TRILINOS_LIBRARY_${_library} CACHE)
  ENDIF()

  FIND_LIBRARY(TRILINOS_LIBRARY_${_library}
    NAMES ${_library}
    HINTS ${Trilinos_LIBRARY_DIRS}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    NO_CMAKE_FIND_ROOT_PATH
    )

  MARK_AS_ADVANCED(TRILINOS_LIBRARY_${_library})

  IF(TRILINOS_LIBRARY_${_library} MATCHES "-NOTFOUND")
    SET(TRILINOS_CONFIG_FOUND FALSE)
  ELSE()
    LIST(APPEND TRILINOS_LIBRARIES ${TRILINOS_LIBRARY_${_library}})
  ENDIF()
ENDFOREACH()

#
# Add the link interface:
#
LIST(APPEND TRILINOS_LIBRARIES
  ${Trilinos_TPL_LIBRARIES}
  ${MPI_CXX_LIBRARIES} # for good measure
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(TRILINOS DEFAULT_MSG
  TRILINOS_LIBRARIES # cosmetic: Gives nice output
  TRILINOS_CONFIG_FOUND
  )

MARK_AS_ADVANCED(TRILINOS_CONFIG EPETRA_CONFIG_H)


IF(TRILINOS_FOUND)
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

  STRING(REGEX REPLACE
    "^[0-9]+\\.[0-9]+\\.([0-9]+).*$" "\\1"
    TRILINOS_VERSION_SUBMINOR "${Trilinos_VERSION}")

  #
  # Determine whether Trilinos was configured with MPI and 64bit indices:
  #
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_MPI_STRING
    REGEX "#define HAVE_MPI")
  IF("${EPETRA_MPI_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_MPI FALSE)
  ELSE()
    SET(TRILINOS_WITH_MPI TRUE)
  ENDIF()
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_32BIT_STRING
    REGEX "#define EPETRA_NO_32BIT_GLOBAL_INDICES")
  IF("${EPETRA_64BIT_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_NO_32BITS_INDICES TRUE)
  ELSE()
    SET(TRILINOS_WITH_NO_32BITS_INDICES FALSE)
  ENDIF()
  FILE(STRINGS "${EPETRA_CONFIG_H}" EPETRA_64BIT_STRING
    REGEX "#define EPETRA_NO_64BIT_GLOBAL_INDICES")
  IF("${EPETRA_64BIT_STRING}" STREQUAL "")
    SET(TRILINOS_WITH_NO_64BITS_INDICES TRUE)
  ELSE()
    SET(TRILINOS_WITH_NO_64BITS_INDICES FALSE)
  ENDIF()

  #
  # Some versions of Sacado_cmath.hpp do things that aren't compatible
  # with the -std=c++0x flag of GCC, see deal.II FAQ.
  # Test whether that is indeed the case:
  #
  IF(_new_trilinos_config)
    UNSET(TRILINOS_SUPPORTS_CPP11 CACHE)
    UNSET(TRILINOS_HAS_C99_TR1_WORKAROUND CACHE)
  ENDIF()

  LIST(APPEND CMAKE_REQUIRED_INCLUDES ${TRILINOS_INCLUDE_DIRS})
  PUSH_TEST_FLAG("${DEAL_II_CXX11_FLAG}")

  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <Sacado_cmath.hpp>
    int main(){ return 0; }
    "
    TRILINOS_SUPPORTS_CPP11
    )

  #
  # Try whether exporting HAS_C99_TR1_CMATH helps:
  #
  PUSH_TEST_FLAG("-DHAS_C99_TR1_CMATH")
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <Sacado_cmath.hpp>
    int main(){ return 0; }
    "
    TRILINOS_HAS_C99_TR1_WORKAROUND
    )

  RESET_CMAKE_REQUIRED()


  MARK_AS_ADVANCED(TRILINOS_DIR)

ELSE()

  SET(TRILINOS_DIR "" CACHE PATH
    "An optional hint to a Trilinos installation"
    )
ENDIF()

