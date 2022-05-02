## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2021 by the deal.II authors
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
#   GMSH_LIBRARY
#   GMSH_INCLUDE_DIR
#   GMSH_WITH_API
#

SET(GMSH_DIR "" CACHE PATH "An optional hint to a Gmsh installation containing the gmsh executable")
SET_IF_EMPTY(GMSH_DIR "$ENV{GMSH_DIR}")

SET(GMSH_LIBRARY_DIR "" CACHE PATH "An optional hint to a Gmsh SDK installation")
SET_IF_EMPTY(GMSH_LIBRARY_DIR "${GMSH_DIR}")

DEAL_II_FIND_PROGRAM(GMSH_EXE gmsh${CMAKE_EXECUTABLE_SUFFIX}
  HINTS ${GMSH_DIR}
  PATH_SUFFIXES bin
  )

DEAL_II_FIND_LIBRARY(GMSH_LIBRARY
  NAMES gmsh
  HINTS ${GMSH_LIBRARY_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_PATH(GMSH_INCLUDE_DIR gmsh.h
  HINTS ${GMSH_LIBRARY_DIR}
  PATH_SUFFIXES include
  )

IF(GMSH_LIBRARY MATCHES "-NOTFOUND" OR 
   GMSH_INCLUDE_DIR MATCHES "-NOTFOUND")
  SET(GMSH_WITH_API FALSE)
ELSE()
  SET(GMSH_WITH_API TRUE)
ENDIF()

DEAL_II_PACKAGE_HANDLE(GMSH
  EXECUTABLE 
    REQUIRED GMSH_EXE
  LIBRARIES
    OPTIONAL GMSH_LIBRARY
  INCLUDE_DIRS 
    OPTIONAL GMSH_INCLUDE_DIR
  USER_INCLUDE_DIRS 
    OPTIONAL GMSH_INCLUDE_DIR
  CLEAR
    GMSH_EXE 
    GMSH_LIBRARY 
    GMSH_INCLUDE_DIR 
    GMSH_WITH_API
  )
