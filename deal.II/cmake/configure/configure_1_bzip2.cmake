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
# Configuration for the bzip2 library:
#

MACRO(FEATURE_BZIP2_FIND_EXTERNAL var)
  FIND_PACKAGE(BZip2)

  IF(BZIP2_FOUND)
    #
    # Rename variables:
    #
    SET(BZIP2_VERSION ${BZIP2_VERSION_STRING})
    SET(BZIP2_INCLUDE_DIRS ${BZIP2_INCLUDE_DIR})

    SET(${var} TRUE)
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(BZIP2)
