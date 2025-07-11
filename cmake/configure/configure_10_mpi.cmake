## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2024 by the deal.II authors
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
  option(DEAL_II_MPI_WITH_DEVICE_SUPPORT "Enable MPI Device support" OFF)
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
