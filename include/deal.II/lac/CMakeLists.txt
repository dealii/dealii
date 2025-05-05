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
    affine_constraints.h
    affine_constraints.templates.h
    arpack_solver.h
    blas_extension_templates.h
    block_indices.h
    block_linear_operator.h
    block_matrix_base.h
    block_sparse_matrix_ez.h
    block_sparse_matrix_ez.templates.h
    block_sparse_matrix.h
    block_sparse_matrix.templates.h
    block_sparsity_pattern.h
    block_vector_base.h
    block_vector.h
    block_vector.templates.h
    chunk_sparse_matrix.h
    chunk_sparse_matrix.templates.h
    chunk_sparsity_pattern.h
    constrained_linear_operator.h
    constraint_matrix.h
    diagonal_matrix.h
    dynamic_sparsity_pattern.h
    eigen.h
    exceptions.h
    full_matrix.h
    full_matrix.templates.h
    generic_linear_algebra.h
    ginkgo_solver.h
    householder.h
    identity_matrix.h
    lapack_full_matrix.h
    lapack_support.h
    lapack_templates.h
    la_parallel_block_vector.h
    la_parallel_block_vector.templates.h
    la_parallel_vector.h
    la_parallel_vector.templates.h
    linear_operator.h
    linear_operator_tools.h
    matrix_block.h
    matrix_iterator.h
    matrix_out.h
    orthogonalization.h
    packaged_operation.h
    parpack_solver.h
    petsc_block_sparse_matrix.h
    petsc_block_vector.h
    petsc_communication_pattern.h
    petsc_compatibility.h
    petsc_full_matrix.h
    petsc_matrix_base.h
    petsc_matrix_free.h
    petsc_parallel_block_sparse_matrix.h
    petsc_parallel_block_vector.h
    petsc_parallel_sparse_matrix.h
    petsc_parallel_vector.h
    petsc_precondition.h
    petsc_snes.h
    petsc_snes.templates.h
    petsc_solver.h
    petsc_sparse_matrix.h
    petsc_ts.h
    petsc_ts.templates.h
    petsc_vector_base.h
    petsc_vector.h
    precondition_block_base.h
    precondition_block.h
    precondition_block.templates.h
    precondition.h
    precondition_selector.h
    qr.h
    read_vector.h
    read_write_vector.h
    read_write_vector.templates.h
    relaxation_block.h
    relaxation_block.templates.h
    scalapack.h
    scalapack.templates.h
    schur_complement.h
    slepc_solver.h
    slepc_spectral_transformation.h
    solver_bicgstab.h
    solver_cg.h
    solver_control.h
    solver_fire.h
    solver_gmres.h
    solver.h
    solver_idr.h
    solver_minres.h
    solver_qmrs.h
    solver_relaxation.h
    solver_richardson.h
    solver_selector.h
    sparse_decomposition.h
    sparse_decomposition.templates.h
    sparse_direct.h
    sparse_ilu.h
    sparse_ilu.templates.h
    sparse_matrix_ez.h
    sparse_matrix_ez.templates.h
    sparse_matrix.h
    sparse_matrix.templates.h
    sparse_matrix_tools.h
    sparse_mic.h
    sparse_mic.templates.h
    sparse_vanka.h
    sparse_vanka.templates.h
    sparsity_pattern_base.h
    sparsity_pattern.h
    sparsity_tools.h
    tensor_product_matrix.h
    tensor_product_matrix.templates.h
    tridiagonal_matrix.h
    trilinos_block_sparse_matrix.h
    trilinos_epetra_communication_pattern.h
    trilinos_epetra_vector.h
    trilinos_index_access.h
    trilinos_linear_operator.h
    trilinos_parallel_block_vector.h
    trilinos_precondition.h
    trilinos_solver.h
    trilinos_sparse_matrix.h
    trilinos_sparsity_pattern.h
    trilinos_tpetra_block_sparse_matrix.h
    trilinos_tpetra_block_sparse_matrix.templates.h
    trilinos_tpetra_block_vector.h
    trilinos_tpetra_block_vector.templates.h
    trilinos_tpetra_communication_pattern.h
    trilinos_tpetra_precondition.h
    trilinos_tpetra_precondition.templates.h
    trilinos_tpetra_solver_direct.h
    trilinos_tpetra_solver_direct.templates.h
    trilinos_tpetra_sparse_matrix.h
    trilinos_tpetra_sparse_matrix.templates.h
    trilinos_tpetra_sparsity_pattern.h
    trilinos_tpetra_types.h
    trilinos_tpetra_vector.h
    trilinos_tpetra_vector.templates.h
    trilinos_vector.h
    utilities.h
    vector_element_access.h
    vector.h
    vector_memory.h
    vector_memory.templates.h
    vector_operation.h
    vector_operations_internal.h
    vector_space_vector.h
    vector.templates.h
    vector_type_traits.h
)


# Prepend each header name by the current directory
list(TRANSFORM _header_files PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Add all of these header files to the global list and propagate things back up one directory level
list(APPEND DEAL_II_HEADER_FILES ${_header_files})
set(DEAL_II_HEADER_FILES ${DEAL_II_HEADER_FILES} PARENT_SCOPE)
