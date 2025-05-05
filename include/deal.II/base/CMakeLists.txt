## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# Configure config.h and revision.h
#
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
  )
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/revision.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/revision.h
  )


#
# Recurse into the sub-directories:
#
add_subdirectory(std_cxx17)
add_subdirectory(std_cxx20)


#
# Keep track of all of the header files in this directory
#
set(_header_files
    aligned_vector.h
    array_view.h
    auto_derivative_function.h
    bounding_box_data_out.h
    bounding_box.h
    communication_pattern_base.h
    complex_overloads.h
    conditional_ostream.h
    convergence_table.h
    data_out_base.h
    derivative_form.h
    discrete_time.h
    enable_observer_pointer.h
    event.h
    exception_macros.h
    exceptions.h
    floating_point_comparator.h
    flow_function.h
    function_bessel.h
    function_cspline.h
    function_derivative.h
    function.h
    function_lib.h
    function_parser.h
    function_restriction.h
    function_signed_distance.h
    function_spherical.h
    function.templates.h
    function_time.h
    function_time.templates.h
    function_tools.h
    geometric_utilities.h
    geometry_info.h
    graph_coloring.h
    hdf5.h
    incremental_function.h
    index_set.h
    init_finalize.h
    iterator_range.h
    job_identifier.h
    kokkos.h
    lazy.h
    linear_index_iterator.h
    logstream.h
    memory_consumption.h
    memory_space_data.h
    memory_space.h
    mg_level_object.h
    mpi_consensus_algorithms.h
    mpi_consensus_algorithms.templates.h
    mpi.h
    mpi_large_count.h
    mpi_noncontiguous_partitioner.h
    mpi_noncontiguous_partitioner.templates.h
    mpi_remote_point_evaluation.h
    mpi_stub.h
    mpi_tags.h
    mpi.templates.h
    multithread_info.h
    mu_parser_internal.h
    mutable_bind.h
    mutex.h
    ndarray.h
    numbers.h
    observer_pointer.h
    parallel.h
    parameter_acceptor.h
    parameter_handler.h
    parsed_convergence_table.h
    parsed_function.h
    partitioner.h
    partitioner.templates.h
    path_search.h
    patterns.h
    point.h
    polynomial.h
    polynomials_abf.h
    polynomials_adini.h
    polynomials_barycentric.h
    polynomials_bdm.h
    polynomials_bernardi_raugel.h
    polynomials_bernstein.h
    polynomials_hermite.h
    polynomials_integrated_legendre_sz.h
    polynomials_nedelec.h
    polynomial_space.h
    polynomials_p.h
    polynomials_piecewise.h
    polynomials_pyramid.h
    polynomials_rannacher_turek.h
    polynomials_raviart_thomas.h
    polynomials_rt_bubbles.h
    polynomials_wedge.h
    process_grid.h
    qprojector.h
    quadrature.h
    quadrature_lib.h
    quadrature_point_data.h
    quadrature_selector.h
    sacado_product_type.h
    scalar_polynomials_base.h
    scope_exit.h
    signaling_nan.h
    smartpointer.h
    subscriptor.h
    symbolic_function.h
    symbolic_function.templates.h
    symmetric_tensor.h
    symmetric_tensor.templates.h
    synchronous_iterator.h
    table.h
    table_handler.h
    table_indices.h
    task_result.h
    template_constraints.h
    tensor_accessors.h
    tensor_function.h
    tensor_function_parser.h
    tensor_function.templates.h
    tensor.h
    tensor_polynomials_base.h
    tensor_product_polynomials_bubbles.h
    tensor_product_polynomials_const.h
    tensor_product_polynomials.h
    thread_local_storage.h
    thread_management.h
    timer.h
    time_stepping.h
    time_stepping.templates.h
    trilinos_utilities.h
    types.h
    undefine_macros.h
    utilities.h
    vectorization.h
    vector_slice.h
    work_stream.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
