#####
##
## Copyright (C) 2012 by the deal.II authors
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
#PROJECT(parameter_gui)

ADD_EXECUTABLE(mesh_conversion
  Main.cc
  MeshConversion.cc
  )

IF("${DEAL_II_EXECUTABLE_RELDIR}" STREQUAL "")
  SET(DEAL_II_EXECUTABLE_RELDIR "bin")
ENDIF()

INSTALL(TARGETS mesh_conversion
  RUNTIME
  DESTINATION ${DEAL_II_EXECUTABLE_RELDIR}
  COMPONENT contrib
  )

