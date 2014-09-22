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
  grid_generator.cc
  grid_in.cc
  grid_out.cc
  grid_refinement.cc
  grid_reordering.cc
  grid_tools.cc
  intergrid_map.cc
  manifold.cc
  manifold_lib.cc
  persistent_tria.cc
  tria_accessor.cc
  tria_boundary.cc
  tria_boundary_lib.cc
  tria.cc
  tria_faces.cc
  tria_levels.cc
  tria_objects.cc
  )

SET(_inst
  grid_generator.inst.in
  grid_in.inst.in
  grid_out.inst.in
  grid_refinement.inst.in
  grid_tools.inst.in
  intergrid_map.inst.in
  manifold.inst.in
  manifold_lib.inst.in
  tria_accessor.inst.in
  tria_boundary.inst.in
  tria_boundary_lib.inst.in
  tria.inst.in
  tria_objects.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/grid/*.h
  )

DEAL_II_ADD_LIBRARY(obj_grid OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_grid "${_inst}")
