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

ADD_LIBRARY(obj_distributed OBJECT ${src})
macro_expand_instantiations(obj_distributed "${inst_in_files}")
