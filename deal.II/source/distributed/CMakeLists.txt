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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

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

