## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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

#
# All source/base files are sufficiently small that they may be added to the
# unity build file (see the developer documentation on DEAL_II_UNITY_BUILD
# for more information).
#
set(_unity_include_src
  auto_derivative_function.cc
  bounding_box.cc
  conditional_ostream.cc
  convergence_table.cc
  discrete_time.cc
  event.cc
  exceptions.cc
  flow_function.cc
  function.cc
  function_cspline.cc
  function_derivative.cc
  function_signed_distance.cc
  function_lib.cc
  function_lib_cutoff.cc
  function_parser.cc
  function_restriction.cc
  function_spherical.cc
  function_time.cc
  function_tools.cc
  geometry_info.cc
  geometric_utilities.cc
  graph_coloring.cc
  incremental_function.cc
  index_set.cc
  job_identifier.cc
  logstream.cc
  hdf5.cc
  kokkos.cc
  mpi.cc
  mpi_compute_index_owner_internal.cc
  mpi_noncontiguous_partitioner.cc
  mpi_remote_point_evaluation.cc
  mu_parser_internal.cc
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
  polynomials_barycentric.cc
  polynomials_bernstein.cc
  polynomials_bdm.cc
  polynomials_bernardi_raugel.cc
  polynomials_hermite.cc
  polynomials_nedelec.cc
  polynomials_integrated_legendre_sz.cc
  polynomial_space.cc
  polynomials_p.cc
  polynomials_piecewise.cc
  polynomials_pyramid.cc
  polynomials_rannacher_turek.cc
  polynomials_raviart_thomas.cc
  polynomials_rt_bubbles.cc
  polynomials_wedge.cc
  process_grid.cc
  qprojector.cc
  quadrature.cc
  quadrature_lib.cc
  quadrature_selector.cc
  scalar_polynomials_base.cc
  subscriptor.cc
  symbolic_function.cc
  table_handler.cc
  tensor.cc
  tensor_function.cc
  tensor_function_parser.cc
  tensor_product_polynomials.cc
  tensor_product_polynomials_bubbles.cc
  tensor_product_polynomials_const.cc
  thread_management.cc
  timer.cc
  time_stepping.cc
  trilinos_utilities.cc
  utilities.cc
  vectorization.cc
  )

set(_separate_src
  data_out_base.cc
  symmetric_tensor.cc
  )

# Add CUDA wrapper files
if(DEAL_II_WITH_CUDA)
  set(_separate_src
    ${_separate_src}
    cuda.cc
  )
endif()

# determined by profiling
set(_n_includes_per_unity_file 29)

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

set(_inst
  bounding_box.inst.in
  data_out_base.inst.in
  function.inst.in
  function_signed_distance.inst.in
  function_restriction.inst.in
  function_time.inst.in
  function_tools.inst.in
  incremental_function.inst.in
  geometric_utilities.inst.in
  mpi.inst.in
  mpi_noncontiguous_partitioner.inst.in
  mpi_remote_point_evaluation.inst.in
  partitioner.inst.in
  polynomials_rannacher_turek.inst.in
  symmetric_tensor.inst.in
  symbolic_function.inst.in
  tensor_function.inst.in
  tensor_function_parser.inst.in
  mu_parser_internal.inst.in
  time_stepping.inst.in
  )

file(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/base/*.h
  )

define_object_library(object_base OBJECT ${_src} ${_header} ${_inst}
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  )
expand_instantiations(object_base "${_inst}")
