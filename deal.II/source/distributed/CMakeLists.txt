SET(src
  grid_refinement.cc
  solution_transfer.cc
  tria.cc
  )

SET(inst_in_files
  grid_refinement.inst.in
  solution_transfer.inst.in
  tria.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_distributed OBJECT ${src})
EXPAND_INSTANTIATIONS(obj_distributed "${inst_in_files}")
DEAL_II_ADD_DEPENDENCIES(obj_distributed deal_ii_target_dependencies)
