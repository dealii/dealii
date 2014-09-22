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
# Configuration for the MUMPS library:
#

SET(FEATURE_MUMPS_DEPENDS MPI LAPACK)


MACRO(FEATURE_MUMPS_FIND_EXTERNAL var)
  FIND_PACKAGE(MUMPS)

  IF(MUMPS_FOUND)
    SET(${var} TRUE)
    CHECK_MPI_INTERFACE(MUMPS ${var})
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(MUMPS)
