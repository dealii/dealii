SET(src
  data_out.cc
  data_out_faces.cc
  data_out_rotation.cc
  data_out_stack.cc
  data_postprocessor.cc
  derivative_approximation.cc
  dof_output_operator.cc
  error_estimator.cc
  fe_field_function.cc
  histogram.cc
  matrix_free.cc
  matrix_tools.cc
  mesh_worker.cc
  mesh_worker_info.cc
  mesh_worker_vector_selector.cc
  point_value_history.cc
  solution_transfer.cc
  time_dependent.cc
  vector_tools.cc
  )

SET(inst_in_files
  data_out_faces.inst.in
  data_out.inst.in
  data_out_rotation.inst.in
  data_out_stack.inst.in
  data_postprocessor.inst.in
  derivative_approximation.inst.in
  dof_output_operator.inst.in
  error_estimator.inst.in
  fe_field_function.inst.in
  matrix_free.inst.in
  matrix_tools.inst.in
  mesh_worker_info.inst.in
  mesh_worker_vector_selector.inst.in
  point_value_history.inst.in
  solution_transfer.inst.in
  time_dependent.inst.in
  vector_tools.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_numerics OBJECT ${src})
EXPAND_INSTANTIATIONS(obj_numerics "${inst_in_files}")
DEAL_II_ADD_DEPENDENCIES(obj_numerics deal_ii_target_dependencies)
