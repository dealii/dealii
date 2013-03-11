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

SET(src
  block_matrix_array.cc
  block_sparse_matrix.cc
  block_sparse_matrix_ez.cc
  block_sparsity_pattern.cc
  block_vector.cc
  chunk_sparse_matrix.cc
  chunk_sparsity_pattern.cc
  compressed_set_sparsity_pattern.cc
  compressed_simple_sparsity_pattern.cc
  compressed_sparsity_pattern.cc
  constraint_matrix.cc
  full_matrix.cc
  lapack_full_matrix.cc
  matrix_lib.cc
  matrix_out.cc
  operator.cc
  parallel_vector.cc
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
  precondition_block.cc
  precondition_block_ez.cc
  relaxation_block.cc
  slepc_solver.cc
  slepc_spectral_transformation.cc
  solver.cc
  solver_control.cc
  sparse_decomposition.cc
  sparse_direct.cc
  sparse_ilu.cc
  sparse_matrix.cc
  sparse_matrix_ez.cc
  sparse_mic.cc
  sparse_vanka.cc
  sparsity_pattern.cc
  sparsity_tools.cc
  swappable_vector.cc
  timestep_control.cc
  tridiagonal_matrix.cc
  trilinos_block_sparse_matrix.cc
  trilinos_block_vector.cc
  trilinos_precondition.cc
  trilinos_solver.cc
  trilinos_sparse_matrix.cc
  trilinos_sparsity_pattern.cc
  trilinos_vector_base.cc
  trilinos_vector.cc
  vector.cc
  vector_memory.cc
  vector_view.cc
)

SET(inst_in_files
  block_sparse_matrix.inst.in
  block_vector.inst.in
  chunk_sparse_matrix.inst.in
  constraint_matrix.inst.in
  full_matrix.inst.in
  lapack_full_matrix.inst.in
  operator.inst.in
  parallel_vector.inst.in
  precondition_block.inst.in
  relaxation_block.inst.in
  sparse_matrix_ez.inst.in
  sparse_matrix.inst.in
  trilinos_sparse_matrix.inst.in
  trilinos_vector_base.inst.in
  vector.inst.in
  vector_memory.inst.in
  vector_view.inst.in
)

DEAL_II_ADD_LIBRARY(obj_lac OBJECT ${src})
EXPAND_INSTANTIATIONS(obj_lac "${inst_in_files}")

