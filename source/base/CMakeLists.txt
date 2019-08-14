## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2019 by the deal.II authors
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

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

#
# All source/base files are sufficiently small that they may be added to the
# unity build file (see the developer documentation on DEAL_II_UNITY_BUILD
# for more information).
#
SET(_unity_include_src
  auto_derivative_function.cc
  bounding_box.cc
  conditional_ostream.cc
  convergence_table.cc
  event.cc
  exceptions.cc
  flow_function.cc
  function.cc
  function_cspline.cc
  function_derivative.cc
  function_lib.cc
  function_lib_cutoff.cc
  function_parser.cc
  function_spherical.cc
  function_time.cc
  geometry_info.cc
  geometric_utilities.cc
  graph_coloring.cc
  incremental_function.cc
  index_set.cc
  job_identifier.cc
  logstream.cc
  hdf5.cc
  mpi.cc
  multithread_info.cc
  named_selection.cc
  numbers.cc
  parallel.cc
  parameter_handler.cc
  parameter_acceptor.cc
  parsed_convergence_table.cc
  parsed_function.cc
  partitioner.cc
  patterns.cc
  path_search.cc
  polynomial.cc
  polynomials_abf.cc
  polynomials_adini.cc
  tensor_polynomials_base.cc
  polynomials_bernstein.cc
  polynomials_bdm.cc
  polynomials_bernardi_raugel.cc
  polynomials_nedelec.cc
  polynomials_integrated_legendre_sz.cc
  polynomial_space.cc
  polynomials_p.cc
  polynomials_piecewise.cc
  polynomials_rannacher_turek.cc
  polynomials_raviart_thomas.cc
  polynomials_rt_bubbles.cc
  process_grid.cc
  quadrature.cc
  quadrature_lib.cc
  quadrature_selector.cc
  subscriptor.cc
  table_handler.cc
  tensor_function.cc
  tensor_product_polynomials.cc
  tensor_product_polynomials_bubbles.cc
  tensor_product_polynomials_const.cc
  thread_management.cc
  tensor.cc
  timer.cc
  time_stepping.cc
  utilities.cc
  vectorization.cc
  )

SET(_separate_src
  data_out_base.cc
  symmetric_tensor.cc
  )

# Add CUDA wrapper files
IF(DEAL_II_WITH_CUDA)
  SET(_separate_src
    ${_separate_src}
    cuda.cu
    partitioner.cu
  )
ENDIF()

# determined by profiling
SET(_n_includes_per_unity_file 29)

SETUP_SOURCE_LIST("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

SET(_inst
  bounding_box.inst.in
  data_out_base.inst.in
  function.inst.in
  function_time.inst.in
  incremental_function.inst.in
  geometric_utilities.inst.in
  hdf5.inst.in
  mpi.inst.in
  partitioner.inst.in
  partitioner.cuda.inst.in
  polynomials_rannacher_turek.inst.in
  symmetric_tensor.inst.in
  tensor.inst.in
  tensor_function.inst.in
  time_stepping.inst.in
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/base/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/base/*.cuh
  )

DEAL_II_ADD_LIBRARY(obj_base OBJECT ${_src} ${_header} ${_inst}
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  )
EXPAND_INSTANTIATIONS(obj_base "${_inst}")
