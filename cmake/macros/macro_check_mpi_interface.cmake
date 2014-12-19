## ---------------------------------------------------------------------
##
## Copyright (C) 2014 by the deal.II authors
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
# Check whether a feature is compiled against the same MPI library as the
# one deal.II picked up
#
# Usage:
#     CHECK_MPI_INTERFACE(_feature _var),
#

MACRO(CHECK_MPI_INTERFACE _feature _var)
  IF(DEAL_II_WITH_MPI)

    SET(_nope FALSE)

    FOREACH(_library ${${_feature}_LIBRARIES})
      IF( _library MATCHES "/libmpi(|_cxx)\\.(a|so)[^/]*$")

        GET_FILENAME_COMPONENT(_file1 ${_library} REALPATH)

        SET(_not_found TRUE)
        FOREACH(_mpi_library ${MPI_LIBRARIES})
          GET_FILENAME_COMPONENT(_file2 ${_mpi_library} REALPATH)
          IF("${_file1}" STREQUAL "${_file2}")
            SET(_not_found FALSE)
            BREAK()
          ENDIF()
        ENDFOREACH()

        IF(_not_found)
          SET(_nope TRUE)
          SET(_spurious_library ${_library})
          BREAK()
        ENDIF()
      ENDIF()
    ENDFOREACH()

    IF(_nope)
      MESSAGE(STATUS "Could not find a sufficient ${_feature} installation: "
        "${_feature} is compiled against a different MPI library than the one "
        "deal.II picked up."
        )
      TO_STRING(_str ${MPI_LIBRARIES})
      SET(PETSC_ADDITIONAL_ERROR_STRING
        ${PETSC_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ${_feature} installation:\n"
        "${_feature} has to be compiled against the same MPI library as deal.II "
        "but the link line of ${_feature} contains:\n"
        "  ${_spurious_library}\n"
        "which is not listed in MPI_LIBRARIES:\n"
        "  MPI_LIBRARIES = \"${_str}\"\n"
        )
      SET(${_var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()

