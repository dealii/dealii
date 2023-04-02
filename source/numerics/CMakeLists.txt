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


# determined by profiling
set(_n_includes_per_unity_file 8)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  cell_data_transfer.inst.in
  data_out_dof_data.inst.in
  data_out_dof_data_codim.inst.in
  data_out_faces.inst.in
  data_out.inst.in
  data_out_resample.inst.in
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
  smoothness_estimator.inst.in
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
  vector_tools_project_codim.inst.in
  vector_tools_project_hp.inst.in
  vector_tools_project_qp.inst.in
  vector_tools_project_qpmf.inst.in
  vector_tools_rhs.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/numerics/*.h
  )

define_object_library(object_numerics OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_numerics "${_inst}")
