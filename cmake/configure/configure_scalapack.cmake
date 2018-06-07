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
# Configuration for the SCALAPACK library:
#

SET(FEATURE_SCALAPACK_DEPENDS MPI LAPACK)


MACRO(FEATURE_SCALAPACK_FIND_EXTERNAL var)
  FIND_PACKAGE(SCALAPACK)

  IF(SCALAPACK_FOUND)
    SET(${var} TRUE)
    CHECK_MPI_INTERFACE(SCALAPACK ${var})
  ENDIF()
ENDMACRO()

CONFIGURE_FEATURE(SCALAPACK)
