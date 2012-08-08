set(SRC_GRID
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

add_library(OBJ_GRID OBJECT ${SRC_GRID})
