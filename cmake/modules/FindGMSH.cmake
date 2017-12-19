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
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find the GMSH library
#
# This module exports
#
#   GMSH_EXECUTABLE
#

SET(GMSH_DIR "" CACHE PATH "An optional hint to a GMSH installation containing the gmsh executable")
SET_IF_EMPTY(GMSH_DIR "$ENV{GMSH_DIR}")

DEAL_II_FIND_FILE(GMSH_EXE gmsh${CMAKE_EXECUTABLE_SUFFIX}
  HINTS ${GMSH_DIR}
  PATH_SUFFIXES bin
  NO_CMAKE_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  )
# NO_CMAKE_PATH and NO_CMAKE_ENVIRONMENT_PATH prevent from looking at
#  <prefix>/include for each <prefix> in CMAKE_PREFIX_PATH before even trying HINTS

DEAL_II_PACKAGE_HANDLE(GMSH
  EXECUTABLE REQUIRED GMSH_EXE
  CLEAR
    GMSH_EXE
  )
