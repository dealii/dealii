## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
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
# Try to find the SCALAPACK library
#
# This module exports
#
#   SCALAPACK_LIBRARIES
#   SCALAPACK_LINKER_FLAGS
#

set(SCALAPACK_DIR "" CACHE PATH "An optional hint to a SCALAPACK directory")
set(BLACS_DIR "" CACHE PATH "An optional hint to a BLACS directory")
set_if_empty(SCALAPACK_DIR "$ENV{SCALAPACK_DIR}")
set_if_empty(BLACS_DIR "$ENV{BLACS_DIR}")

#
# Search for scalapack:
#

deal_ii_find_library(SCALAPACK_LIBRARY NAMES scalapack scalapack-openmpi scalapack-pvm scalapack-mpi scalapack-mpich scalapack-mpich2 scalapack-lam
  HINTS ${SCALAPACK_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Well, depending on the version of scalapack and the distribution it might
# be necessary to search for blacs, too.
if (SCALAPACK_LIBRARY)
  message(STATUS "Check if BLACS is embedded in ScaLAPACK library")

  clear_cmake_required()
  set(CMAKE_REQUIRED_FLAGS ${DEAL_II_CXX_FLAGS_SAVED} ${MPI_CXX_COMPILE_FLAGS} ${MPI_CXX_LINK_FLAGS})
  set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_PATH})
  set(CMAKE_REQUIRED_LIBRARIES ${DEAL_II_LINKER_FLAGS_SAVED} ${MPI_LIBRARIES} ${SCALAPACK_LIBRARY} ${LAPACK_LIBRARIES})
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
  reset_cmake_required()

  # If Blacs is not embedded, try to find it
  if (NOT SCALAPACK_LIBRARY_HAS_BLACS)
    message(STATUS "Try to find BLACS")
    # FIXME: use instead:
    #   find_package(BLACS ${SCALAPACK_FIND_QUIETLY})
    #   set(SCALAPACK_LIBRARY ${SCALAPACK_LIBRARY} ${BLACS_LIBRARY})

    # FIXME: this is potentially dangerous endeavour as we can easily pickup
    # system-provided BLACS build with some MPI which is inconsistent with
    # ScaLAPACK library found above.
    foreach(_lib blacs blacsCinit blacsF77init)
      string(TOUPPER "${_lib}" _lib_upper)
      deal_ii_find_library(${_lib_upper}_LIBRARY
        NAMES ${_lib} ${_lib}_MPI-LINUX-0 ${_lib}_MPI-DARWIN-0 ${_lib}-openmpi
        HINTS ${BLACS_DIR} ${SCALAPACK_DIR} ${SCALAPACK_DIR}/../blacs/
        PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib LIB
      )
    endforeach()
  endif()
endif()


process_feature(SCALAPACK
  LIBRARIES
    REQUIRED SCALAPACK_LIBRARY LAPACK_LIBRARIES
    OPTIONAL BLACS_LIBRARY BLACSCINIT_LIBRARY BLACSF77INIT_LIBRARY MPI_Fortran_LIBRARIES
  LINKER_FLAGS
    OPTIONAL LAPACK_LINKER_FLAGS
  CLEAR
    SCALAPACK_LIBRARY
    BLACS_LIBRARY BLACSCINIT_LIBRARY BLACSF77INIT_LIBRARY
  )
