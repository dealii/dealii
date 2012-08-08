set(SRC_NUMERICS
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
  operator.cc
  point_value_history.cc
  solution_transfer.cc
  time_dependent.cc
  vectors.cc
  )

add_library(OBJ_NUMERICS OBJECT ${SRC_NUMERICS})
