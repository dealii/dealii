## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
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
  evaluation_template_factory.cc
  evaluation_template_factory_inst2.cc
  evaluation_template_factory_inst3.cc
  evaluation_template_factory_inst4.cc
  evaluation_template_factory_inst5.cc
  evaluation_template_factory_inst6.cc
  evaluation_template_factory_hanging_nodes.cc
  evaluation_template_face_factory.cc
  evaluation_template_face_factory_inst2.cc
  evaluation_template_face_factory_inst3.cc
  evaluation_template_face_factory_inst4.cc
  evaluation_template_face_factory_inst5.cc
  evaluation_template_face_factory_inst6.cc
  fe_point_evaluation.cc
  mapping_info.cc
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

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/matrix_free/*.h
  )

define_object_library(object_matrix_free OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_matrix_free "${_inst}")


# Keep track of a complete list of all source files, based on the
# variables we have set above.  Prepend each source name by the
# current directory.
set(_source_files)
list(APPEND _source_files ${_src})
list(TRANSFORM _source_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these source files to the global list and propagate
# things back up one directory level
list(APPEND DEAL_II_SOURCE_FILES ${_source_files})
set(DEAL_II_SOURCE_FILES ${DEAL_II_SOURCE_FILES} PARENT_SCOPE)
