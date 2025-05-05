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
# Keep track of all of the header files in this directory
#
set(_header_files
    adaptation_strategies.h
    cell_data_transfer.h
    cell_data_transfer.templates.h
    data_component_interpretation.h
    data_out_dof_data.h
    data_out_dof_data.templates.h
    data_out_faces.h
    data_out.h
    data_out_resample.h
    data_out_rotation.h
    data_out_stack.h
    data_postprocessor.h
    derivative_approximation.h
    dof_output_operator.h
    dof_output_operator.templates.h
    dof_print_solver_step.h
    error_estimator.h
    error_estimator.templates.h
    fe_field_function.h
    fe_field_function.templates.h
    histogram.h
    history.h
    matrix_creator.h
    matrix_creator.templates.h
    matrix_tools.h
    nonlinear.h
    point_value_history.h
    rtree.h
    smoothness_estimator.h
    solution_transfer.h
    tensor_product_matrix_creator.h
    time_dependent.h
    vector_tools_boundary.h
    vector_tools_boundary.templates.h
    vector_tools_common.h
    vector_tools_constraints.h
    vector_tools_constraints.templates.h
    vector_tools_evaluate.h
    vector_tools.h
    vector_tools_integrate_difference.h
    vector_tools_integrate_difference.templates.h
    vector_tools_interpolate.h
    vector_tools_interpolate.templates.h
    vector_tools_mean_value.h
    vector_tools_mean_value.templates.h
    vector_tools_point_gradient.h
    vector_tools_point_gradient.templates.h
    vector_tools_point_value.h
    vector_tools_point_value.templates.h
    vector_tools_project.h
    vector_tools_project.templates.h
    vector_tools_rhs.h
    vector_tools_rhs.templates.h
    vector_tools.templates.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
