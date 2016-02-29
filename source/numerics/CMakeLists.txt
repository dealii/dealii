## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  data_out.cc
  data_out_dof_data.cc
  data_out_faces.cc
  data_out_rotation.cc
  data_out_stack.cc
  data_postprocessor.cc
  derivative_approximation.cc
  dof_output_operator.cc
  error_estimator_1d.cc
  error_estimator.cc
  error_estimator_inst2.cc
  fe_field_function.cc
  histogram.cc
  matrix_creator.cc
  matrix_creator_inst2.cc
  matrix_creator_inst3.cc
  matrix_tools_once.cc
  matrix_tools.cc
  point_value_history.cc
  solution_transfer.cc
  solution_transfer_inst2.cc
  solution_transfer_inst3.cc
  solution_transfer_inst4.cc
  time_dependent.cc
  vector_tools_boundary.cc
  vector_tools_constraints.cc
  vector_tools_integrate_difference.cc
  vector_tools_interpolate.cc
  vector_tools_mean_value.cc
  vector_tools_point_value.cc
  vector_tools_point_gradient.cc
  vector_tools_project.cc
  vector_tools_rhs.cc
  )

SET(_inst
  data_out_dof_data.inst.in
  data_out_faces.inst.in
  data_out.inst.in
  data_out_rotation.inst.in
  data_out_stack.inst.in
  data_postprocessor.inst.in
  derivative_approximation.inst.in
  dof_output_operator.inst.in
  error_estimator_1d.inst.in
  error_estimator.inst.in
  fe_field_function.inst.in
  matrix_creator.inst.in
  matrix_tools.inst.in
  point_value_history.inst.in
  solution_transfer.inst.in
  time_dependent.inst.in
  vector_tools_boundary.inst.in
  vector_tools_constraints.inst.in
  vector_tools_integrate_difference.inst.in
  vector_tools_interpolate.inst.in
  vector_tools_mean_value.inst.in
  vector_tools_point_value.inst.in
  vector_tools_point_gradient.inst.in
  vector_tools_project.inst.in
  vector_tools_rhs.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/numerics/*.h
  )

DEAL_II_ADD_LIBRARY(obj_numerics OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_numerics "${_inst}")
