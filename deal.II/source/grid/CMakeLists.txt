#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

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

DEAL_II_ADD_LIBRARY(obj_grid OBJECT ${src})
EXPAND_INSTANTIATIONS(obj_grid "${inst_in_files}")

