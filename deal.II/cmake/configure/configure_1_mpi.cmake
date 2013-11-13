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
# Configuration for mpi support:
#
# We look for the C and Fortran libraries as well because they are needed
# by some external libraries:
#

MACRO(FEATURE_MPI_FIND_EXTERNAL var)
  #
  # Obey a manual user override: If MPI_CXX_FOUND is set to true in the
  # cache, we skip the FIND_PACKAGE calls:
  #
  IF(MPI_CXX_FOUND)
    SET(MPI_FOUND TRUE)
  ENDIF()

  #
  # If CMAKE_CXX_COMPILER is already an MPI wrapper, use it to determine
  # the mpi implementation. If MPI_CXX_COMPILER is defined use the value
  # directly.
  #
  SET_IF_EMPTY(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
  IF(CMAKE_C_COMPILER_WORKS)
    SET_IF_EMPTY(MPI_C_COMPILER ${CMAKE_C_COMPILER}) # for good measure
  ELSE()
    MESSAGE(STATUS
      "No suitable C compiler was found! MPI C interface can not be "
      "autodetected"
      )
  ENDIF()
  IF(CMAKE_Fortran_COMPILER_WORKS)
    SET_IF_EMPTY(MPI_Fortran_COMPILER ${CMAKE_Fortran_COMPILER}) # for good measure
  ELSE()
    MESSAGE(STATUS
      "No suitable Fortran compiler was found! MPI Fortran interface can "
      "not be autodetected"
      )
  ENDIF()
  FIND_PACKAGE(MPI)

  IF(NOT MPI_CXX_FOUND AND DEAL_II_WITH_MPI)
    #
    # CMAKE_CXX_COMPILER is apparently not an mpi wrapper.
    # So, let's be a bit more aggressive in finding MPI (and if
    # DEAL_II_WITH_MPI is set).
    #
    MESSAGE(STATUS
      "MPI not found but DEAL_II_WITH_MPI is set to TRUE."
      " Try again with more aggressive search paths:"
      )
    SET(MPI_FOUND) # clear this value so that FIND_PACKAGE runs again.
    UNSET(MPI_CXX_COMPILER CACHE)
    UNSET(MPI_C_COMPILER CACHE)
    UNSET(MPI_Fortran_COMPILER CACHE)
    FIND_PACKAGE(MPI)
  ENDIF()

  #
  # Manually clean up variables:
  #
  FOREACH(_lang C CXX Fortran)
    IF(MPI_${_lang}_LIBRARIES MATCHES "-NOTFOUND")
      SET(MPI_${_lang}_LIBRARIES)
    ENDIF()
  ENDFOREACH()

  IF(MPI_CXX_FOUND)
    #
    # Manually assemble some version information:
    #
    FIND_FILE(MPI_MPI_H NAMES mpi.h
      HINTS ${MPI_INCLUDE_PATH}
      )

    IF(NOT MPI_MPI_H MATCHES "-NOTFOUND" AND NOT DEFINED MPI_VERSION)
      FILE(STRINGS "${MPI_MPI_H}" MPI_VERSION_MAJOR_STRING
        REGEX "#define.*MPI_VERSION")
      STRING(REGEX REPLACE "^.*MPI_VERSION.*([0-9]+).*" "\\1"
        MPI_VERSION_MAJOR "${MPI_VERSION_MAJOR_STRING}"
        )
      FILE(STRINGS ${MPI_MPI_H} MPI_VERSION_MINOR_STRING
        REGEX "#define.*MPI_SUBVERSION")
      STRING(REGEX REPLACE "^.*MPI_SUBVERSION.*([0-9]+).*" "\\1"
        MPI_VERSION_MINOR "${MPI_VERSION_MINOR_STRING}"
        )
      SET(MPI_VERSION "${MPI_VERSION_MAJOR}.${MPI_VERSION_MINOR}")
      IF("${MPI_VERSION}" STREQUAL ".")
        SET(MPI_VERSION)
        SET(MPI_VERSION_MAJOR)
        SET(MPI_VERSION_MINOR)
      ENDIF()

      # OMPI specific version number:
      FILE(STRINGS ${MPI_MPI_H} OMPI_VERSION_MAJOR_STRING
        REGEX "#define.*OMPI_MAJOR_VERSION")
      STRING(REGEX REPLACE "^.*OMPI_MAJOR_VERSION.*([0-9]+).*" "\\1"
        OMPI_VERSION_MAJOR "${OMPI_VERSION_MAJOR_STRING}"
        )
      FILE(STRINGS ${MPI_MPI_H} OMPI_VERSION_MINOR_STRING
        REGEX "#define.*OMPI_MINOR_VERSION")
      STRING(REGEX REPLACE "^.*OMPI_MINOR_VERSION.*([0-9]+).*" "\\1"
        OMPI_VERSION_MINOR "${OMPI_VERSION_MINOR_STRING}"
        )
      FILE(STRINGS ${MPI_MPI_H} OMPI_VERSION_RELEASE_STRING
        REGEX "#define.*OMPI_RELEASE_VERSION")
      STRING(REGEX REPLACE "^.*OMPI_RELEASE_VERSION.*([0-9]+).*" "\\1"
        OMPI_VERSION_SUBMINOR "${OMPI_VERSION_RELEASE_STRING}"
        )
      SET(OMPI_VERSION
        "${OMPI_VERSION_MAJOR}.${OMPI_VERSION_MINOR}.${OMPI_VERSION_SUBMINOR}"
        )
      IF("${OMPI_VERSION}" STREQUAL "..")
        SET(OMPI_VERSION)
        SET(OMPI_VERSION_MAJOR)
        SET(OMPI_VERSION_MINOR)
        SET(OMPI_VERSION_SUBMINOR)
      ENDIF()
    ENDIF()

    SET(${var} TRUE)
  ENDIF()

  # Hide some variables:
  MARK_AS_ADVANCED(MPI_EXTRA_LIBRARY MPI_LIBRARY MPI_MPI_H)

ENDMACRO()


MACRO(FEATURE_MPI_CONFIGURE_EXTERNAL)
  #
  # The user has to know the location of the mpi headers as well:
  #
  SET(MPI_CXX_ADD_TO_USER_INCLUDE_DIRS TRUE)

  REGISTER_FEATURE(MPI_CXX)
ENDMACRO()


MACRO(FEATURE_MPI_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable mpi library!\n"
    "Please ensure that an mpi library is installed on your computer\n"
    "and set CMAKE_CXX_COMPILER to the appropriate mpi wrappers:\n"
    "    $ CXX=\".../mpicxx\" cmake <...>\n"
    "    $ cmake -DCMAKE_CXX_COMPILER=\".../mpicxx\" <...>\n"
    "Or with additional C and Fortran wrappers (recommended!):\n"
    "    $ CC=\".../mpicc\" CXX=\".../mpicxx\" F90=\".../mpif90\" cmake <...>\n"
    "    $ cmake -DCMAKE_C_COMPILER=\".../mpicc\"\\\n"
    "            -DCMAKE_CXX_COMPILER=\".../mpicxx\"\\\n"
    "            -DCMAKE_Fortran_COMPILER=\".../mpif90\"\\\n"
    "            <...>\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(MPI)
