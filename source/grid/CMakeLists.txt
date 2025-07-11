## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Work around an issue where gcc emits a note that the warning
# -Wmisleading-indentation will be disabled due to a humongous line of over
# 100k characters produced by one of the CGAL headers. Unfortunately,
# guarding the include by DEAL_II_DISABLE_EXTRA_DIAGNOSTICS is not enough
# due to a longstanding bug in gcc. Thus, simply set
# -Wno-misleading-indentation on the command line for the
# grid_generator_cgal.cc and grid_tools.cc compilation units.
# See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89549
#
if(DEAL_II_WITH_CGAL)
  enable_if_supported(_flag "-Wno-misleading-indentation")
  set_property(SOURCE "grid_generator_cgal.cc" "grid_tools.cc"
    APPEND PROPERTY COMPILE_OPTIONS "${_flag}"
    )
endif()

set(_unity_include_src
  cell_id.cc
  grid_refinement.cc
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

set(_separate_src
  grid_generator.cc
  grid_generator_cgal.cc
  grid_generator_from_name.cc
  grid_generator_pipe_junction.cc
  grid_in.cc
  grid_out.cc
  grid_tools.cc
  grid_tools_cache.cc
  grid_tools_dof_handlers.cc
  grid_tools_geometry.cc
  grid_tools_nontemplates.cc
  grid_tools_topology.cc
  tria.cc
  )

# include all files in the unity file
set(_n_includes_per_unity_file 20)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  cell_id.inst.in
  grid_generator_cgal.inst.in
  grid_generator_from_name.inst.in
  grid_generator.inst.in
  grid_generator_pipe_junction.inst.in
  grid_in.inst.in
  grid_out.inst.in
  grid_refinement.inst.in
  grid_tools.inst.in
  grid_tools_cache.inst.in
  grid_tools_dof_handlers.inst.in
  grid_tools_geometry.inst.in
  grid_tools_topology.inst.in
  intergrid_map.inst.in
  manifold.inst.in
  manifold_lib.inst.in
  reference_cell.inst.in
  tria_accessor.inst.in
  tria_description.inst.in
  tria.inst.in
  tria_objects.inst.in
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/grid/*.h
  )

define_object_library(object_grid OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_grid "${_inst}")

collect_source_and_header_files("${_unity_include_src};${_separate_src}" "${_header}")
