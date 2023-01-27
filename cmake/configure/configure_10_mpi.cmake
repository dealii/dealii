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
set(DEAL_II_WITH_MPI OFF CACHE BOOL "")

#
# Configuration for mpi support:
#

macro(feature_mpi_find_external var)
  find_package(DEAL_II_MPI)

  if(MPI_FOUND)
    set(${var} TRUE)

    if(NOT MPI_HAVE_MPI_SEEK_SET)
      message(STATUS
        "Could not find a sufficient MPI version: "
        "Your MPI implementation must define MPI_SEEK_SET.")
      set(MPI_ADDITIONAL_ERROR_STRING
        "Your MPI implementation must define MPI_SEEK_SET.\n")
      set(${var} FALSE)
    endif()

    if(MPI_VERSION VERSION_LESS "3.0")
      message(STATUS
        "Could not find a sufficient MPI version: "
        "Your MPI implementation does not support the MPI 3.0 standard.")
      set(MPI_ADDITIONAL_ERROR_STRING
        "Your MPI implementation does not support the MPI 3.0 standard.\n")
      set(${var} FALSE)
    endif()

  endif()
endmacro()

macro(feature_mpi_configure_external)

  #
  # We must convert the MPIEXEC_(PRE|POST)FLAGS strings to lists in order
  # to use them in command lines:
  #
  separate_arguments(MPIEXEC_PREFLAGS)
  separate_arguments(MPIEXEC_POSTFLAGS)

  #
  # TODO: We might consider refactoring this option into an automatic check
  # (in Modules/FindMPI.cmake) at some point. For the time being this is an
  # advanced configuration option.
  #
  if(DEAL_II_MPI_WITH_CUDA_SUPPORT)
    option(DEAL_II_MPI_WITH_DEVICE_SUPPORT "Enable MPI Device support" ON)
  else()
    option(DEAL_II_MPI_WITH_DEVICE_SUPPORT "Enable MPI Device support" OFF)
  endif()
  mark_as_advanced(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
endmacro()

macro(feature_mpi_error_message)
  message(FATAL_ERROR "\n"
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
endmacro()


configure_feature(MPI)


if(NOT DEAL_II_WITH_MPI)
  #
  # Disable and hide the DEAL_II_MPI_WITH_DEVICE_SUPPORT option
  #
  set(DEAL_II_MPI_WITH_DEVICE_SUPPORT)
  unset(DEAL_II_MPI_WITH_DEVICE_SUPPORT CACHE)
endif()
