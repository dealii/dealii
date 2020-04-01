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

SET(_unity_include_src
  affine_constraints.cc
  block_sparse_matrix.cc
  block_sparse_matrix_ez.cc
  block_sparsity_pattern.cc
  block_vector.cc
  chunk_sparse_matrix.cc
  chunk_sparsity_pattern.cc
  dynamic_sparsity_pattern.cc
  exceptions.cc
  scalapack.cc
  la_vector.cc
  la_parallel_vector.cc
  la_parallel_block_vector.cc
  matrix_out.cc
  precondition_block.cc
  precondition_block_ez.cc
  relaxation_block.cc
  read_write_vector.cc
  solver.cc
  solver_bicgstab.cc
  solver_control.cc
  sparse_decomposition.cc
  sparse_direct.cc
  sparse_ilu.cc
  sparse_matrix_ez.cc
  sparse_mic.cc
  sparse_vanka.cc
  sparsity_pattern.cc
  sparsity_tools.cc
  vector.cc
  vector_memory.cc
  )

# make sure that files requiring the LAPACK wrappers are compiled separately.
SET(_separate_src
  full_matrix.cc
  lapack_full_matrix.cc
  sparse_matrix.cc
  sparse_matrix_inst2.cc
  tridiagonal_matrix.cc
  )

# Add CUDA wrapper files
IF(DEAL_II_WITH_CUDA)
  SET(_separate_src
    ${_separate_src}
    vector_memory.cu
    affine_constraints.cu
    )
ENDIF()

SET(_inst
  affine_constraints.inst.in
  block_sparse_matrix.inst.in
  block_vector.inst.in
  chunk_sparse_matrix.inst.in
  full_matrix.inst.in
  lapack_full_matrix.inst.in
  la_vector.inst.in
  la_parallel_vector.inst.in
  la_parallel_block_vector.inst.in
  precondition_block.inst.in
  relaxation_block.inst.in
  read_write_vector.inst.in
  scalapack.inst.in
  solver.inst.in
  sparse_matrix_ez.inst.in
  sparse_matrix.inst.in
  vector.inst.in
  vector_memory.inst.in
  vector_memory_release.inst.in
  )


# Add PETSc wrapper files. If PETSc has not been found,
# then these files should be empty and there is no need
# to even look at them
IF(DEAL_II_WITH_PETSC)
  SET(_unity_include_src
    ${_unity_include_src}
    petsc_full_matrix.cc
    petsc_matrix_base.cc
    petsc_matrix_free.cc
    petsc_parallel_block_sparse_matrix.cc
    petsc_parallel_block_vector.cc
    petsc_parallel_sparse_matrix.cc
    petsc_parallel_vector.cc
    petsc_precondition.cc
    petsc_solver.cc
    petsc_sparse_matrix.cc
    petsc_vector_base.cc
  )
ENDIF()

# Same for SLEPc
IF(DEAL_II_WITH_SLEPC)
  SET(_unity_include_src
    ${_unity_include_src}
    slepc_solver.cc
    slepc_spectral_transformation.cc
  )
ENDIF()

IF(DEAL_II_WITH_GINKGO)
  SET(_unity_include_src
    ${_unity_include_src}
    ginkgo_solver.cc
  )
ENDIF()

# Also add Trilinos wrapper files
IF(DEAL_II_WITH_TRILINOS)
  SET(_unity_include_src
    ${_unity_include_src}
    trilinos_block_sparse_matrix.cc
    trilinos_block_vector.cc
    trilinos_epetra_communication_pattern.cc
    trilinos_epetra_vector.cc
    trilinos_precondition.cc
    trilinos_precondition_ml.cc
    trilinos_solver.cc
    trilinos_sparse_matrix.cc
    trilinos_sparsity_pattern.cc
    trilinos_tpetra_communication_pattern.cc
    trilinos_tpetra_vector.cc
    trilinos_vector.cc
  )

  SET(_separate_src
    ${_separate_src}
    trilinos_precondition_muelu.cc
    )

  SET(_inst
    ${_inst}
    trilinos_sparse_matrix.inst.in
    trilinos_vector.inst.in
  )
ENDIF()

# Add CUDA wrapper files
IF(DEAL_II_WITH_CUDA)
  SET(_separate_src
    ${_separate_src}
    cuda_kernels.cu
    cuda_precondition.cu
    cuda_solver_direct.cu
    cuda_sparse_matrix.cu
    cuda_vector.cu
    la_parallel_vector.cu
    read_write_vector.cu
  )
ENDIF()

# determined by profiling
SET(_n_includes_per_unity_file 20)

IF(DEAL_II_UNITY_BUILD)
  # sort files so that the petsc and trilinos files are not together
  LIST(SORT _unity_include_src)
ENDIF()

SETUP_SOURCE_LIST("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/lac/*.h
  )

DEAL_II_ADD_LIBRARY(obj_lac OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_lac "${_inst}")
