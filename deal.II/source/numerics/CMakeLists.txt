#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/numerics/*.h
  )

SET(_src
  data_out.cc
  data_out_faces.cc
  data_out_rotation.cc
  data_out_stack.cc
  data_postprocessor.cc
  derivative_approximation.cc
  dof_output_operator.cc
  error_estimator.cc
  fe_field_function.cc
  matrix_tools.cc
  histogram.cc
  point_value_history.cc
  solution_transfer.cc
  solution_transfer_inst2.cc
  time_dependent.cc
  vector_tools.cc
  )

SET(_inst
  data_out_faces.inst.in
  data_out.inst.in
  data_out_rotation.inst.in
  data_out_stack.inst.in
  data_postprocessor.inst.in
  derivative_approximation.inst.in
  dof_output_operator.inst.in
  error_estimator.inst.in
  fe_field_function.inst.in
  matrix_tools.inst.in
  point_value_history.inst.in
  solution_transfer.inst.in
  time_dependent.inst.in
  vector_tools.inst.in
  )

DEAL_II_ADD_LIBRARY(obj_numerics OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_numerics "${_inst}")

