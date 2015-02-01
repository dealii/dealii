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
# Try to find the ARPACK library
#
# This module exports
#
#   ARPACK_LIBRARIES
#   ARPACK_LINKER_FLAGS
#

#
# TODO: ARPACK and mpi...
#

SET(ARPACK_DIR "" CACHE PATH "An optional hint to an ARPACK installation")
SET_IF_EMPTY(ARPACK_DIR "$ENV{ARPACK_DIR}")

DEAL_II_FIND_LIBRARY(ARPACK_LIBRARY
  NAMES arpack
  HINTS ${ARPACK_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_PACKAGE_HANDLE(ARPACK
  LIBRARIES REQUIRED ARPACK_LIBRARY LAPACK_LIBRARIES
  LINKER_FLAGS OPTIONAL LAPACK_LINKER_FLAGS
  CLEAR ARPACK_LIBRARY
  )
