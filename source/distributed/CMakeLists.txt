## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set(_unity_include_src
  grid_refinement.cc
  cell_weights.cc
  cell_data_transfer.cc
  fully_distributed_tria.cc
  repartitioning_policy_tools.cc
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
  tria.inst.in
  shared_tria.inst.in
  tria_base.inst.in
  p4est_wrappers.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/distributed/*.h
  )

define_object_library(object_distributed OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_distributed "${_inst}")
