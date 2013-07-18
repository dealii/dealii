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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
#PROJECT(mesh_converter)

ADD_EXECUTABLE(mesh_converter
  Main.cc
  MeshConversion.cc
  )

IF("${DEAL_II_EXECUTABLE_RELDIR}" STREQUAL "")
  SET(DEAL_II_EXECUTABLE_RELDIR "bin")
ENDIF()

INSTALL(TARGETS mesh_converter
  RUNTIME
  DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
  COMPONENT mesh_converter
  )

