## ---------------------------------------------------------------------
## $Id$
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
      IF(_library MATCHES "/libmpi[^/]*\\.so")
        LIST(FIND MPI_LIBRARIES "${_library}" _position)
        IF("${_position}" STREQUAL "-1")
          SET(_nope TRUE)
          SET(_mpi_library ${_library})
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
        "  ${_mpi_library}\n"
        "which is not listed in MPI_LIBRARIES:\n"
        "  MPI_LIBRARIES = \"${_str}\"\n"
        )
      SET(${_var} FALSE)
    ENDIF()
  ENDIF()
ENDMACRO()

