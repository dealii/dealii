## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# Configuration for the p4est and sc libraries:
#

SET(FEATURE_P4EST_DEPENDS MPI)


MACRO(FEATURE_P4EST_FIND_EXTERNAL var)
  FIND_PACKAGE(P4EST)

  IF(P4EST_FOUND)
    SET(${var} TRUE)

    #
    # We require at least version 0.3.4.1
    #
    IF(P4EST_VERSION VERSION_LESS  "0.3.4.1")
      MESSAGE(STATUS "Insufficient p4est installation found: "
        "At least version 0.3.4.1 is required."
        )
      SET(P4EST_ADDITIONAL_ERROR_STRING
        "Insufficient p4est installation found!\n"
        "At least version 0.3.4.1 is required.\n"
        )
      SET(${var} FALSE)
    ENDIF()

    #
    # Check whether p4est supports mpi:
    #
    IF(NOT P4EST_WITH_MPI)
      MESSAGE(STATUS "Insufficient p4est installation found: "
        "p4est has to be configured with MPI enabled."
        )
      SET(P4EST_ADDITIONAL_ERROR_STRING
        ${P4EST_ADDITIONAL_ERROR_STRING}
        "Insufficient p4est installation found!\n"
        "p4est has to be configured with MPI enabled.\n"
        )
      SET(${var} FALSE)
    ENDIF()

    CHECK_MPI_INTERFACE(P4EST ${var})
  ENDIF()
ENDMACRO()


CONFIGURE_FEATURE(P4EST)
