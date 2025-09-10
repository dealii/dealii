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
  mg_base.cc
  mg_constrained_dofs.cc
  mg_level_global_transfer.cc
  mg_transfer_block.cc
  mg_transfer_component.cc
  mg_transfer_internal.cc
  mg_transfer_prebuilt.cc
  multigrid.cc
  )

set(_separate_src
  mg_tools.cc
  mg_transfer_global_coarsening.cc
  mg_transfer_matrix_free.cc
  )

# concatenate all unity inclusion files in one file
set(_n_includes_per_unity_file 15)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  mg_base.inst.in
  mg_constrained_dofs.inst.in
  mg_level_global_transfer.inst.in
  mg_tools.inst.in
  mg_transfer_block.inst.in
  mg_transfer_component.inst.in
  mg_transfer_global_coarsening.inst.in
  mg_transfer_internal.inst.in
  mg_transfer_matrix_free.inst.in
  mg_transfer_prebuilt.inst.in
  multigrid.inst.in
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/multigrid/*.h
  )

define_object_library(object_multigrid OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_multigrid "${_inst}")

collect_source_and_header_files("${_unity_include_src};${_separate_src}" "${_header}")
