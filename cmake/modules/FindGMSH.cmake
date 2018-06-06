## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2018 by the deal.II authors
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
# Try to find the GMSH library
#
# This module exports
#
#   GMSH_EXECUTABLE
#

SET(GMSH_DIR "" CACHE PATH "An optional hint to a Gmsh installation containing the gmsh executable")
SET_IF_EMPTY(GMSH_DIR "$ENV{GMSH_DIR}")

DEAL_II_FIND_PROGRAM(GMSH_EXE gmsh${CMAKE_EXECUTABLE_SUFFIX}
  HINTS ${GMSH_DIR}
  PATH_SUFFIXES bin
  )

DEAL_II_PACKAGE_HANDLE(GMSH
  EXECUTABLE REQUIRED GMSH_EXE
  CLEAR
    GMSH_EXE
  )
