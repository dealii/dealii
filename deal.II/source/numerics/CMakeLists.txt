SET(src
  data_out.cc
  data_out_faces.cc
  data_out_rotation.cc
  data_out_stack.cc
  data_postprocessor.cc
  derivative_approximation.cc
  error_estimator.cc
  fe_field_function.cc
  histogram.cc
  matrices.cc
  matrix_free.cc
  mesh_worker.cc
  mesh_worker_info.cc
  mesh_worker_vector_selector.cc
  dof_output_operator.cc
  point_value_history.cc
  solution_transfer.cc
  time_dependent.cc
  vectors.cc
  )

SET(inst_in_files
  data_out_faces.inst.in
  data_out.inst.in
  data_out_rotation.inst.in
  data_out_stack.inst.in
  data_postprocessor.inst.in
  derivative_approximation.inst.in
  error_estimator.inst.in
  fe_field_function.inst.in
  matrices.inst.in
  matrix_free.inst.in
  mesh_worker_info.inst.in
  mesh_worker_vector_selector.inst.in
  dof_output_operator.inst.in
  point_value_history.inst.in
  solution_transfer.inst.in
  time_dependent.inst.in
  vectors.inst.in
  )

ADD_LIBRARY(obj_numerics OBJECT ${src})
macro_expand_instantiations(obj_numerics "${inst_in_files}")
ADD_DEPENDENCIES(obj_numerics deal_ii_target_dependencies)
