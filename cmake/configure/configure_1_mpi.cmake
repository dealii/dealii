## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2016 by the deal.II authors
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

    CHECK_COMPILER_SETUP(
      "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_SAVED} ${MPI_CXX_FLAGS}"
      "${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_SAVED} ${MPI_LINKER_FLAGS}"
      MPI_WORKING_COMPILER
      ${DEAL_II_LIBRARIES} ${MPI_LIBRARIES}
      )

    IF(NOT MPI_WORKING_COMPILER)
      #
      # Try a workaround and drop "-fuse-ld=gold" (if present) from the
      # linker invocation
      #
      MESSAGE(STATUS "Unable to compile a simple test program. "
        "Try to drop \"-fuse-ld=gold\" from the linker flags."
        )
      STRING(REPLACE "-fuse-ld=gold" "" _filtered_flags "${DEAL_II_LINKER_FLAGS}")

      CHECK_COMPILER_SETUP(
        "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_SAVED} ${MPI_CXX_FLAGS}"
        "${_filtered_flags} ${DEAL_II_LINKER_FLAGS_SAVED} ${MPI_LINKER_FLAGS}"
        MPI_WORKING_COMPILER
        ${DEAL_II_LIBRARIES} ${MPI_LIBRARIES}
        )

      IF(MPI_WORKING_COMPILER)
        SET(DEAL_II_LINKER_FLAGS "${_filtered_flags}")
      ELSE()
        MESSAGE(STATUS "Could not find a sufficient MPI installation: "
          "Unable to compile a simple test program."
          )
        SET(MPI_ADDITIONAL_ERROR_STRING
          ${MPI_ADDITIONAL_ERROR_STRING}
          "Unable to compile and link a simple test program with your MPI installation. \n"
          )
        SET(${var} FALSE)
      ENDIF()
    ENDIF()

  ENDIF()
ENDMACRO()

MACRO(FEATURE_MPI_ERROR_MESSAGE)
  MESSAGE(FATAL_ERROR "\n"
    "Could not find any suitable mpi library!\n"
    ${MPI_ADDITIONAL_ERROR_STRING}
    "\nPlease ensure that an mpi library is installed on your computer\n"
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
