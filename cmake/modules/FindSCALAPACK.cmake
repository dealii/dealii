## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
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
# Try to find the SCALAPACK library
#
# This module exports
#
#   SCALAPACK_LIBRARIES
#   SCALAPACK_LINKER_FLAGS
#

SET(SCALAPACK_DIR "" CACHE PATH "An optional hint to a SCALAPACK directory")
SET(BLACS_DIR "" CACHE PATH "An optional hint to a BLACS directory")
SET_IF_EMPTY(SCALAPACK_DIR "$ENV{SCALAPACK_DIR}")
SET_IF_EMPTY(BLACS_DIR "$ENV{BLACS_DIR}")

#
# Search for scalapack:
#

DEAL_II_FIND_LIBRARY(SCALAPACK_LIBRARY NAMES scalapack scalapack-openmpi scalapack-pvm scalapack-mpi scalapack-mpich scalapack-mpich2 scalapack-lam
  HINTS ${SCALAPACK_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Well, depending on the version of scalapack and the distribution it might
# be necessary to search for blacs, too.
IF (SCALAPACK_LIBRARY)
  MESSAGE(STATUS "Check if BLACS is embedded in ScaLAPACK library")

  CLEAR_CMAKE_REQUIRED()
  SET(CMAKE_REQUIRED_FLAGS ${MPI_CXX_COMPILE_FLAGS} ${MPI_CXX_LINK_FLAGS})
  SET(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
  SET(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES} ${SCALAPACK_LIBRARY} ${LAPACK_LIBRARIES})
  CHECK_CXX_SOURCE_COMPILES(
    "
    #include <mpi.h>
    extern \"C\"
    {
      int Csys2blacs_handle(MPI_Comm comm);
    }
    int main() {
      const int res = Csys2blacs_handle(MPI_COMM_WORLD);
    }
    "
    SCALAPACK_LIBRARY_HAS_BLACS
  )
  RESET_CMAKE_REQUIRED()

  # If Blacs is not embedded, try to find it
  IF (NOT SCALAPACK_LIBRARY_HAS_BLACS)
    MESSAGE(STATUS "Try to find BLACS")
    # FIXME: use instead:
    #   FIND_PACKAGE(BLACS ${SCALAPACK_FIND_QUIETLY})
    #   SET(SCALAPACK_LIBRARY ${SCALAPACK_LIBRARY} ${BLACS_LIBRARY})

    # FIXME: this is potentially dangerous endeavour as we can easily pickup
    # system-provided BLACS build with some MPI which is inconsistent with
    # ScaLAPACK library found above.
    FOREACH(_lib blacs blacsCinit blacsF77init)
      STRING(TOUPPER "${_lib}" _lib_upper)
      DEAL_II_FIND_LIBRARY(${_lib_upper}_LIBRARY
        NAMES ${_lib} ${_lib}_MPI-LINUX-0 ${_lib}_MPI-DARWIN-0 ${_lib}-openmpi
        HINTS ${BLACS_DIR} ${SCALAPACK_DIR} ${SCALAPACK_DIR}/../blacs/
        PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib LIB
      )
    ENDFOREACH()
  ENDIF()
ENDIF()


DEAL_II_PACKAGE_HANDLE(SCALAPACK
  LIBRARIES
    REQUIRED SCALAPACK_LIBRARY LAPACK_LIBRARIES
    OPTIONAL BLACS_LIBRARY BLACSCINIT_LIBRARY BLACSF77INIT_LIBRARY MPI_Fortran_LIBRARIES
  LINKER_FLAGS
    OPTIONAL LAPACK_LINKER_FLAGS
  CLEAR
    SCALAPACK_LIBRARY
    BLACS_LIBRARY BLACSCINIT_LIBRARY BLACSF77INIT_LIBRARY
  )
