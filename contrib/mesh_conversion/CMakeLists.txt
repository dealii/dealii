## ---------------------------------------------------------------------
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

IF("${DEAL_II_EXECUTABLE_RELDIR}" STREQUAL "")
  SET(DEAL_II_EXECUTABLE_RELDIR "bin")
ENDIF()

ADD_EXECUTABLE(mesh_converter_exe
  Main.cc
  MeshConversion.cc
  )
SET_TARGET_PROPERTIES(mesh_converter_exe
  PROPERTIES
  OUTPUT_NAME mesh_converter
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${DEAL_II_EXECUTABLE_RELDIR}"
  )

ADD_DEPENDENCIES(mesh_converter mesh_converter_exe)

INSTALL(TARGETS mesh_converter_exe
  RUNTIME DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
  COMPONENT mesh_converter
  )
EXPORT(TARGETS mesh_converter_exe
  FILE
  ${CMAKE_BINARY_DIR}/${DEAL_II_PROJECT_CONFIG_RELDIR}/${DEAL_II_PROJECT_CONFIG_NAME}Targets.cmake
  APPEND
  )
