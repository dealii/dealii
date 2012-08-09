SET(src
  grid_generator.cc
  grid_in.cc
  grid_out.cc
  grid_refinement.cc
  grid_reordering.cc
  grid_tools.cc
  intergrid_map.cc
  persistent_tria.cc
  tria_accessor.cc
  tria_boundary.cc
  tria_boundary_lib.cc
  tria.cc
  tria_faces.cc
  tria_levels.cc
  tria_objects.cc
  )

SET(inst_in_files
  grid_generator.inst.in
  grid_in.inst.in
  grid_out.inst.in
  grid_refinement.inst.in
  grid_tools.inst.in
  intergrid_map.inst.in
  tria_accessor.inst.in
  tria_boundary.inst.in
  tria_boundary_lib.inst.in
  tria.inst.in
  tria_objects.inst.in
  )

ADD_LIBRARY(obj_grid OBJECT ${src})
macro_expand_instantiations(obj_grid "${inst_in_files}")
ADD_DEPENDENCIES(obj_grid contrib)
