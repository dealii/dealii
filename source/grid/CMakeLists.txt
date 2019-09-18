## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2018 by the deal.II authors
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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_unity_include_src
  cell_id.cc
  grid_generator.cc
  grid_generator_from_name.cc
  grid_in.cc
  grid_out.cc
  grid_refinement.cc
  grid_tools_cache.cc
  intergrid_map.cc
  manifold.cc
  manifold_lib.cc
  persistent_tria.cc
  tria_accessor.cc
  tria_faces.cc
  tria_levels.cc
  tria_objects.cc
  )

SET(_separate_src
  grid_reordering.cc
  grid_tools.cc
  grid_tools_dof_handlers.cc
  tria.cc
  )

# include all files in the unity file
SET(_n_includes_per_unity_file 20)

SETUP_SOURCE_LIST("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

SET(_inst
  cell_id.inst.in
  grid_generator.inst.in
  grid_generator_from_name.inst.in
  grid_in.inst.in
  grid_out.inst.in
  grid_refinement.inst.in
  grid_tools.inst.in
  grid_tools_dof_handlers.inst.in
  grid_tools_cache.inst.in
  intergrid_map.inst.in
  manifold.inst.in
  manifold_lib.inst.in
  tria_accessor.inst.in
  tria.inst.in
  tria_objects.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/grid/*.h
  )

DEAL_II_ADD_LIBRARY(obj_grid OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_grid "${_inst}")
