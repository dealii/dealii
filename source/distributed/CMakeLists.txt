## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2021 by the deal.II authors
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

include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

set(_unity_include_src
  grid_refinement.cc
  cell_weights.cc
  cell_data_transfer.cc
  fully_distributed_tria.cc
  repartitioning_policy_tools.cc
  solution_transfer.cc
  tria.cc
  tria_base.cc
  shared_tria.cc
  p4est_wrappers.cc
  field_transfer.cc
  )

set(_separate_src
  )

# concatenate all unity inclusion files in one file
set(_n_includes_per_unity_file 15)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  grid_refinement.inst.in
  cell_weights.inst.in
  cell_data_transfer.inst.in
  field_transfer.inst.in
  fully_distributed_tria.inst.in
  repartitioning_policy_tools.inst.in
  solution_transfer.inst.in
  tria.inst.in
  shared_tria.inst.in
  tria_base.inst.in
  p4est_wrappers.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/distributed/*.h
  )

deal_ii_add_library(obj_distributed OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(obj_distributed "${_inst}")
