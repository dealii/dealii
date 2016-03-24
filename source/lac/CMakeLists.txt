## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_CURRENT_BINARY_DIR})

SET(_src
  block_matrix_array.cc
  block_sparse_matrix.cc
  block_sparse_matrix_ez.cc
  block_sparsity_pattern.cc
  block_vector.cc
  chunk_sparse_matrix.cc
  chunk_sparsity_pattern.cc
  dynamic_sparsity_pattern.cc
  constraint_matrix.cc
  full_matrix.cc
  lapack_full_matrix.cc
  la_vector.cc
  matrix_lib.cc
  matrix_out.cc
  parallel_vector.cc
  precondition_block.cc
  precondition_block_ez.cc
  relaxation_block.cc
  read_write_vector.cc
  solver.cc
  solver_control.cc
  sparse_decomposition.cc
  sparse_direct.cc
  sparse_ilu.cc
  sparse_matrix.cc
  sparse_matrix_inst2.cc
  sparse_matrix_ez.cc
  sparse_mic.cc
  sparse_vanka.cc
  sparsity_pattern.cc
  sparsity_tools.cc
  swappable_vector.cc
  tridiagonal_matrix.cc
  vector.cc
  vector_memory.cc
  vector_view.cc
  )

SET(_inst
  block_sparse_matrix.inst.in
  block_vector.inst.in
  chunk_sparse_matrix.inst.in
  constraint_matrix.inst.in
  full_matrix.inst.in
  lapack_full_matrix.inst.in
  la_vector.inst.in
  parallel_vector.inst.in
  precondition_block.inst.in
  relaxation_block.inst.in
  read_write_vector.inst.in
  solver.inst.in
  sparse_matrix_ez.inst.in
  sparse_matrix.inst.in
  vector.inst.in
  vector_memory.inst.in
  vector_view.inst.in
  )

# Add PETSc wrapper files. If PETSc has not been found,
# then these files should be empty and there is no need
# to even look at them
IF(DEAL_II_WITH_PETSC)
  SET(_src
    ${_src}
    petsc_block_sparse_matrix.cc
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
    petsc_vector.cc
  )
ENDIF()

# Same for SLEPc
IF(DEAL_II_WITH_SLEPC)
  SET(_src
    ${_src}
    slepc_solver.cc
    slepc_spectral_transformation.cc
  )
ENDIF()

# Also add Trilinos wrapper files
IF(DEAL_II_WITH_TRILINOS)
  SET(_src
    ${_src}
    trilinos_block_sparse_matrix.cc
    trilinos_block_vector.cc
    trilinos_epetra_communication_pattern.cc
    trilinos_epetra_vector.cc
    trilinos_parallel_block_vector.cc
    trilinos_precondition.cc
    trilinos_precondition_ml.cc
    trilinos_precondition_muelu.cc
    trilinos_solver.cc
    trilinos_sparse_matrix.cc
    trilinos_sparsity_pattern.cc
    trilinos_vector_base.cc
    trilinos_vector.cc
  )
  SET(_inst
    ${_inst}
    trilinos_sparse_matrix.inst.in
    trilinos_vector_base.inst.in
  )
ENDIF()



FILE(GLOB _header
  ${CMAKE_SOURCE_DIR}/include/deal.II/lac/*.h
  )

DEAL_II_ADD_LIBRARY(obj_lac OBJECT ${_src} ${_header} ${_inst})
EXPAND_INSTANTIATIONS(obj_lac "${_inst}")
