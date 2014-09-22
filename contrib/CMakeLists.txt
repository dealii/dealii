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

IF(DEAL_II_COMPONENT_MESH_CONVERTER)
  MESSAGE(STATUS "Setting up mesh_converter")
  ADD_SUBDIRECTORY(mesh_conversion)
  MESSAGE(STATUS "Setting up mesh_converter - Done")
ENDIF()

IF(DEAL_II_COMPONENT_PARAMETER_GUI)
  MESSAGE(STATUS "Setting up parameter_gui")
  ADD_SUBDIRECTORY(parameter_gui)
  MESSAGE(STATUS "Setting up parameter_gui - Done")
ENDIF()

