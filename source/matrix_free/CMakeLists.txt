## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set(_src
  dof_info.cc
  evaluation_template_factory_inst1.cc
  evaluation_template_factory_inst2.cc
  evaluation_template_factory_inst3.cc
  evaluation_template_factory_inst4.cc
  evaluation_template_factory_inst5.cc
  evaluation_template_factory_inst6.cc
  evaluation_template_factory_hanging_nodes.cc
  evaluation_template_face_factory_inst1.cc
  evaluation_template_face_factory_inst2.cc
  evaluation_template_face_factory_inst3.cc
  evaluation_template_face_factory_inst4.cc
  evaluation_template_face_factory_inst5.cc
  evaluation_template_face_factory_inst6.cc
  fe_point_evaluation.cc
  mapping_info_inst1.cc
  mapping_info_inst2.cc
  mapping_info_inst3.cc
  matrix_free.cc
  portable_matrix_free.cc
  shape_info.cc
  task_info.cc
  vector_data_exchange.cc
  )

set(_inst
  evaluation_template_factory.inst.in
  evaluation_template_face_factory.inst.in
  evaluation_template_factory_hanging_nodes.inst.in
  fe_point_evaluation.inst.in
  mapping_info.inst.in
  matrix_free.inst.in
  shape_info.inst.in
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/matrix_free/*.h
  )

define_object_library(object_matrix_free OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_matrix_free "${_inst}")

collect_source_and_header_files("${_src}" "${_header}")
