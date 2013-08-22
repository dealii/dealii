## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
# Configuration for the metis library:
#

MACRO(FEATURE_METIS_FIND_EXTERNAL var)
  FIND_PACKAGE(METIS)

  IF(METIS_FOUND)
    IF(METIS_VERSION_MAJOR GREATER 4)
      SET(${var} TRUE)
    ELSE()
      MESSAGE(STATUS "Insufficient metis installation found: "
        "Version 5.x required!"
        )
      SET(METIS_ADDITIONAL_ERROR_STRING
        "Could not find a sufficient modern metis installation: "
        "Version 5.x required!\n"
        )

      UNSET(METIS_LIBRARY CACHE)
      UNSET(METIS_INCLUDE_DIR CACHE)
      SET(METIS_DIR "" CACHE PATH
        "An optional hint to a metis directory"
        )
      MARK_AS_ADVANCED(CLEAR METIS_DIR)
    ENDIF()
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(METIS)
