## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/dofs/*.h
  )

define_object_library(object_dofs OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_dofs "${_inst}")

collect_source_and_header_files("${_unity_include_src};${_separate_src}" "${_header}")
