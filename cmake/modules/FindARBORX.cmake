## ---------------------------------------------------------------------
##
## Copyright (C) 2021 by the deal.II authors
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
# Try to find the ArborX library
#
# This module exports
#
#   ARBORX_INCLUDE_DIRS
#   ARBORX_INTERFACE_LINK_FLAGS
#

SET(ARBORX_DIR "" CACHE PATH "An optional hint to an ArborX installation")
SET_IF_EMPTY(ARBORX_DIR "$ENV{ARBORX_DIR}")


FIND_PACKAGE(ArborX
  HINTS ${ARBORX_DIR} ${ArborX_DIR} $ENV{ArborX_DIR}
  )

SET(_libraries "")

IF(ArborX_FOUND)
  GET_PROPERTY(ARBORX_INSTALL_INCLUDE_DIR TARGET ArborX::ArborX PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
ENDIF()

DEAL_II_PACKAGE_HANDLE(ARBORX
  LIBRARIES REQUIRED ${_libraries}
  INCLUDE_DIRS REQUIRED ARBORX_INSTALL_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED ARBORX_INSTALL_INCLUDE_DIR
  CLEAR ARBORX_DIR ${_libraries}
  )
