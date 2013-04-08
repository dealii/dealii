#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is multi licensed under
##   - QPL 1.0,
##   - GPL 2.0 or any later version,
##   - LGL 2.1 or any later version.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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

