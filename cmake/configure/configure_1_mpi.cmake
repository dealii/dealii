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
# Never autodetect MPI support. This forces the user to explicitly enable
# MPI support with -DWITH_MPI=ON (or -DDEAL_II_WITH_MPI=ON) on the command
# line.
#
SET(DEAL_II_WITH_MPI OFF CACHE BOOL "")

#
# Configuration for mpi support:
#

MACRO(FEATURE_MPI_FIND_EXTERNAL var)
  FIND_PACKAGE(MPI)

  IF(MPI_FOUND)
    SET(${var} TRUE)

    IF(NOT MPI_HAVE_MPI_SEEK_SET)
      MESSAGE(STATUS
        "Could not find a sufficient MPI version: "
        "Your MPI implementation must define MPI_SEEK_SET.")
      SET(MPI_ADDITIONAL_ERROR_STRING
        "Your MPI implementation must define MPI_SEEK_SET.\n")
      SET(${var} FALSE)
    ENDIF()

    IF(MPI_VERSION VERSION_LESS "3.0")
      MESSAGE(STATUS
        "Could not find a sufficient MPI version: "
        "Your MPI implementation does not support the MPI 3.0 standard.")
      SET(MPI_ADDITIONAL_ERROR_STRING
        "Your MPI implementation does not support the MPI 3.0 standard.\n")
      SET(${var} FALSE)
    ENDIF()

  ENDIF()
ENDMACRO()

MACRO(FEATURE_MPI_CONFIGURE_EXTERNAL)

  #
  # We must convert the MPIEXEC_(PRE|POST)FLAGS strings to lists in order
  # to use them in command lines:
  #
  SEPARATE_ARGUMENTS(MPIEXEC_PREFLAGS)
  SEPARATE_ARGUMENTS(MPIEXEC_POSTFLAGS)

  #
  # TODO: We might consider refactoring this option into an automatic check
  # (in Modules/FindMPI.cmake) at some point. For the time being this is an
  # advanced configuration option.
  #
  OPTION(DEAL_II_MPI_WITH_CUDA_SUPPORT "Enable MPI Cuda support" OFF)
  MARK_AS_ADVANCED(DEAL_II_MPI_WITH_CUDA_SUPPORT)
ENDMACRO()

MACRO(FEATURE_MPI_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable mpi library!\n"
    ${MPI_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that an mpi library is installed on your computer\n"
    "and set MPI_CXX_COMPILER to the appropriate mpi wrappers:\n"
    "    $ cmake -DMPI_CXX_COMPILER=\".../mpicxx\" <...>\n"
    "Or with additional C and Fortran wrappers (recommended!):\n"
    "    $ cmake -DMPI_C_COMPILER=\".../mpicc\"\\\n"
    "            -DMPI_CXX_COMPILER=\".../mpicxx\"\\\n"
    "            -DMPI_Fortran_COMPILER=\".../mpif90\"\\\n"
    "            <...>\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(MPI)


IF(NOT DEAL_II_WITH_MPI)
  #
  # Disable and hide the DEAL_II_MPI_WITH_CUDA_SUPPORT option
  #
  SET(DEAL_II_MPI_WITH_CUDA_SUPPORT)
  UNSET(DEAL_II_MPI_WITH_CUDA_SUPPORT CACHE)
ENDIF()
