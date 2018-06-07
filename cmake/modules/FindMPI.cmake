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
# Find MPI
#
# This module exports:
#   MPI_LIBRARIES
#   MPI_INCLUDE_DIRS
#   MPI_CXX_FLAGS
#   MPI_LINKER_FLAGS
#   MPI_VERSION
#   MPI_VERSION_MAJOR
#   MPI_VERSION_MINOR
#   OMPI_VERSION
#   MPI_HAVE_MPI_SEEK_SET
#

#
# Configuration for mpi support:
#
# We look for the C and Fortran libraries as well because they are needed
# by some external libraries for the link interface.
#

IF(MPI_CXX_FOUND)
  SET(MPI_FOUND TRUE)
ENDIF()

#
# Call the system FindMPI.cmake module:
#

# in case MPIEXEC is specified first call find_program() so that in case of
# success its subsequent runs inside FIND_PACKAGE(MPI) do not alter the
# desired result.
IF(DEFINED ENV{MPIEXEC})
  FIND_PROGRAM(MPIEXEC $ENV{MPIEXEC})
ENDIF()

# temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)
FIND_PACKAGE(MPI)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

#
# Older versions of MPI may not have MPI_SEEK_SET, which we
# require. Strangely, unlike MPICH, OpenMPI needs the correct link libraries
# for this to compile, not *just* the correct include directories.
#

CLEAR_CMAKE_REQUIRED()
SET(CMAKE_REQUIRED_FLAGS ${MPI_CXX_COMPILE_FLAGS} ${MPI_CXX_LINK_FLAGS})
SET(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
SET(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <mpi.h>
  #ifndef MPI_SEEK_SET
  #  error
  #endif
  int main() {}
  "
  MPI_HAVE_MPI_SEEK_SET
  )
RESET_CMAKE_REQUIRED()

#
# Newer versions of FindMPI.cmake only populate MPI_CXX_* (and MPI_C_*,
# MPI_Fortran_*) variables. Let's rename these version names
#

IF(NOT DEFINED MPI_VERSION AND DEFINED MPI_CXX_VERSION)
  SET(MPI_VERSION ${MPI_CXX_VERSION})
  SET(MPI_VERSION_MAJOR ${MPI_CXX_VERSION_MAJOR})
  SET(MPI_VERSION_MINOR ${MPI_CXX_VERSION_MINOR})
ENDIF()

#
# Really old versions of CMake do not export any version information. In
# this case, query the mpi.h header for the necessary information:
#

DEAL_II_FIND_FILE(MPI_MPI_H
  NAMES mpi.h
  HINTS ${MPI_CXX_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH}
  )
IF(NOT MPI_MPI_H MATCHES "-NOTFOUND" AND NOT DEFINED MPI_VERSION)
  FILE(STRINGS "${MPI_MPI_H}" MPI_VERSION_MAJOR_STRING
    REGEX "#define.*MPI_VERSION")
  STRING(REGEX REPLACE "^.*MPI_VERSION[ ]+([0-9]+).*" "\\1"
    MPI_VERSION_MAJOR "${MPI_VERSION_MAJOR_STRING}"
    )
  FILE(STRINGS ${MPI_MPI_H} MPI_VERSION_MINOR_STRING
    REGEX "#define.*MPI_SUBVERSION")
  STRING(REGEX REPLACE "^.*MPI_SUBVERSION[ ]+([0-9]+).*" "\\1"
    MPI_VERSION_MINOR "${MPI_VERSION_MINOR_STRING}"
    )
  SET(MPI_VERSION "${MPI_VERSION_MAJOR}.${MPI_VERSION_MINOR}")
ENDIF()

#
# Except - this doesn't always work. Some distributions install a header
# stub mpi.h that includes the right mpi header depending on the
# architecture. In this case we are really out of luck. It is not
# straightforward to find the correct header file to query the version
# information from. Just set a very conservative default:
#
IF(NOT DEFINED MPI_VERSION OR MPI_VERSION STREQUAL ".")
  SET(MPI_VERSION "0.0")
  SET(MPI_VERSION_MAJOR "0")
  SET(MPI_VERSION_MINOR "0")
ENDIF()

DEAL_II_PACKAGE_HANDLE(MPI
  LIBRARIES
    OPTIONAL MPI_CXX_LIBRARIES MPI_Fortran_LIBRARIES MPI_C_LIBRARIES
  INCLUDE_DIRS
    OPTIONAL MPI_CXX_INCLUDE_PATH MPI_C_INCLUDE_PATH
  USER_INCLUDE_DIRS
    OPTIONAL MPI_CXX_INCLUDE_PATH MPI_C_INCLUDE_PATH
  CXX_FLAGS OPTIONAL MPI_CXX_COMPILE_FLAGS
  LINKER_FLAGS OPTIONAL MPI_CXX_LINK_FLAGS
  CLEAR
    MPI_C_COMPILER
    MPI_CXX_COMPILER
    MPIEXEC
    MPI_EXTRA_LIBRARY
    MPI_Fortran_COMPILER
    MPI_HEADER_PATH
    MPI_LIB
    MPI_LIBRARY
    MPI_MPI_H
    MPI_HAVE_MPI_SEEK_SET
  )

