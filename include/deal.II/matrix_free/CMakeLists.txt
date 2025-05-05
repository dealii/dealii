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
    constraint_info.h
    dof_info.h
    dof_info.templates.h
    evaluation_flags.h
    evaluation_kernels_face.h
    evaluation_kernels.h
    evaluation_kernels_common.h
    evaluation_kernels_hanging_nodes.h
    evaluation_selector.h
    evaluation_template_face_factory.templates.h
    evaluation_template_factory.h
    evaluation_template_factory_hanging_nodes.templates.h
    evaluation_template_factory_internal.h
    evaluation_template_factory.templates.h
    face_info.h
    face_setup_internal.h
    fe_evaluation_data.h
    fe_evaluation.h
    fe_point_evaluation.h
    fe_remote_evaluation.h
    hanging_nodes_internal.h
    mapping_data_on_the_fly.h
    mapping_info.h
    mapping_info_storage.h
    mapping_info_storage.templates.h
    mapping_info.templates.h
    matrix_free.h
    matrix_free.templates.h
    operators.h
    portable_evaluation_kernels.h
    portable_fe_evaluation.h
    portable_hanging_nodes_internal.h
    portable_matrix_free.h
    portable_matrix_free.templates.h
    portable_tensor_product_kernels.h
    shape_info.h
    shape_info.templates.h
    task_info.h
    tensor_product_kernels.h
    tensor_product_point_kernels.h
    tools.h
    type_traits.h
    util.h
    vector_access_internal.h
    vector_data_exchange.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
