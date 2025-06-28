## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
  enable_observer_pointer.cc
  event.cc
  exceptions.cc
  flow_function.cc
  function.cc
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
  init_finalize.cc
  index_set.cc
  job_identifier.cc
  logstream.cc
  kokkos.cc
  mpi.cc
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
  polynomials_vector_anisotropic.cc
  polynomials_wedge.cc
  qprojector.cc
  quadrature.cc
  quadrature_lib.cc
  quadrature_selector.cc
  scalar_polynomials_base.cc
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

if(DEAL_II_WITH_SCALAPACK)
  set(_unity_include_src
    ${_unity_include_src}
    process_grid.cc
  )
endif()

if(DEAL_II_WITH_HDF5)
  set(_unity_include_src
    ${_unity_include_src}
    hdf5.cc
  )
endif()

if(DEAL_II_WITH_GSL)
  set(_unity_include_src
    ${_unity_include_src}
    function_cspline.cc
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


# Collect the header files inside the base/ header directory. There
# are some headers in sub-folders as well that are not represented by
# a corresponding directory under source/base/, so glob them here as
# well.
#
# Finally, there are also header files in
# include/deal.II/boost_adaptors/ that would be otherwise unaccounted
# for, and we'll collect them here as well, in the absence of a better
# place.
file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/base/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/base/std_cxx*/*.h
  ${CMAKE_SOURCE_DIR}/include/deal.II/boost_adaptors/*.h
  )

define_object_library(object_base OBJECT ${_src} ${_header} ${_inst}
  ${CMAKE_BINARY_DIR}/include/deal.II/base/config.h
  )
expand_instantiations(object_base "${_inst}")

collect_source_and_header_files("${_unity_include_src};${_separate_src}" "${_header}")
