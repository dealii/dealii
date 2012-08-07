set(SRC_DISTRIBUTED
  grid_refinement.cc
  solution_transfer.cc
  tria.cc
  )

add_library(OBJ_DISTRIBUTED OBJECT ${SRC_DISTRIBUTED})
