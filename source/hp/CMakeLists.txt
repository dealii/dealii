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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  dof_faces.cc
  dof_handler.cc
  dof_level.cc
  fe_collection.cc
  fe_values.cc
  mapping_collection.cc
  )

SET(_inst
  dof_handler.inst.in
  fe_collection.inst.in
  fe_values.inst.in
  mapping_collection.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/hp/*.h
  )

DEAL_II_ADD_LIBRARY(obj_hp OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_hp "${_inst}")
