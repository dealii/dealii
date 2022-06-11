## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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

#
# Work around an issue where gcc emits a note that the warning
# -Wmisleading-indentation will be disabled due to a humongous line of over
# 100k characters produced by one of the CGAL headers. Unfortunately,
# guarding the include by DEAL_II_DISABLE_EXTRA_DIAGNOSTICS is not enough
# due to a longstanding bug in gcc. Thus, simply wet
# -Wno-misleading-indentation on the command line for the
# grid_generator_cgal.cc compilation unit.
#
IF(DEAL_II_WITH_CGAL)
  ENABLE_IF_SUPPORTED(_flag "-Wno-misleading-indentation")
  SET_PROPERTY(SOURCE "grid_generator_cgal.cc"
    APPEND_STRING PROPERTY COMPILE_FLAGS " ${_flag}"
    )
ENDIF()


SET(_unity_include_src
  cell_id.cc
  grid_refinement.cc
  grid_reordering.cc
  intergrid_map.cc
  manifold.cc
  manifold_lib.cc
  persistent_tria.cc
  reference_cell.cc
  tria_accessor.cc
  tria_description.cc
  tria_faces.cc
  tria_levels.cc
  tria_objects.cc
  )

SET(_separate_src
  grid_generator.cc
  grid_generator_cgal.cc
  grid_generator_from_name.cc
  grid_generator_pipe_junction.cc
  grid_in.cc
  grid_out.cc
  grid_tools_cache.cc
  grid_tools.cc
  grid_tools_dof_handlers.cc
  grid_tools_nontemplates.cc
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
  grid_generator_cgal.inst.in
  grid_generator_from_name.inst.in
  grid_generator.inst.in
  grid_generator_pipe_junction.inst.in
  grid_in.inst.in
  grid_out.inst.in
  grid_refinement.inst.in
  grid_tools.inst.in
  grid_tools_dof_handlers.inst.in
  grid_tools_cache.inst.in
  intergrid_map.inst.in
  manifold.inst.in
  manifold_lib.inst.in
  reference_cell.inst.in
  tria_accessor.inst.in
  tria_description.inst.in
  tria.inst.in
  tria_objects.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/grid/*.h
  )

DEAL_II_ADD_LIBRARY(obj_grid OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_grid "${_inst}")
