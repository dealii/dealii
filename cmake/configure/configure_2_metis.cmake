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

SET(FEATURE_METIS_AFTER MPI)

#
# Configuration for the metis library:
#

MACRO(FEATURE_METIS_FIND_EXTERNAL var)
  FIND_PACKAGE(METIS)

  IF(METIS_FOUND)
    SET(${var} TRUE)

    IF(NOT METIS_VERSION_MAJOR GREATER 4)
      MESSAGE(STATUS "Insufficient metis installation found: "
        "Version 5.x required!"
        )
      SET(METIS_ADDITIONAL_ERROR_STRING
        "Could not find a sufficient modern metis installation: "
        "Version 5.x required!\n"
        )
      SET(${var} FALSE)
    ENDIF()

    CHECK_MPI_INTERFACE(METIS ${var})
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(METIS)
