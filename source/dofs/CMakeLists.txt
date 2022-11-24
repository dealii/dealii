## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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
  block_info.cc
  dof_faces.cc
  dof_handler.cc
  dof_objects.cc
  number_cache.cc
  )

set(_separate_src
  dof_accessor.cc
  dof_accessor_get.cc
  dof_accessor_set.cc
  dof_handler_policy.cc
  dof_renumbering.cc
  dof_tools.cc
  dof_tools_constraints.cc
  dof_tools_sparsity.cc
  )

# concatenate all unity inclusion files in one file
set(_n_includes_per_unity_file 15)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  block_info.inst.in
  dof_accessor_get.inst.in
  dof_accessor.inst.in
  dof_accessor_set.inst.in
  dof_handler.inst.in
  dof_handler_policy.inst.in
  dof_objects.inst.in
  dof_renumbering.inst.in
  dof_tools_constraints.inst.in
  dof_tools.inst.in
  dof_tools_sparsity.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/dofs/*.h
  )

deal_ii_add_library(obj_dofs OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(obj_dofs "${_inst}")
