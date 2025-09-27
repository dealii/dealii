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

set(_unity_include_src
  affine_constraints.cc
  block_sparse_matrix.cc
  block_sparse_matrix_ez.cc
  block_sparsity_pattern.cc
  block_vector.cc
  chunk_sparse_matrix.cc
  chunk_sparsity_pattern.cc
  dynamic_sparsity_pattern.cc
  exceptions.cc
  la_parallel_vector.cc
  la_parallel_block_vector.cc
  matrix_out.cc
  precondition_block.cc
  precondition_block_ez.cc
  relaxation_block.cc
  read_write_vector.cc
  solver.cc
  solver_control.cc
  solver_gmres.cc
  sparse_decomposition.cc
  sparse_direct.cc
  sparse_ilu.cc
  sparse_matrix_ez.cc
  sparse_mic.cc
  sparse_vanka.cc
  sparsity_pattern_base.cc
  sparsity_pattern.cc
  sparsity_tools.cc
  tensor_product_matrix.cc
  vector.cc
  vector_memory.cc
  )

# make sure that files requiring the LAPACK wrappers are compiled separately.
set(_separate_src
  full_matrix.cc
  lapack_full_matrix.cc
  qr.cc
  sparse_matrix_inst1.cc
  sparse_matrix_inst2.cc
  tridiagonal_matrix.cc
  utilities.cc
  )

set(_inst
  affine_constraints.inst.in
  block_sparse_matrix.inst.in
  block_vector.inst.in
  chunk_sparse_matrix.inst.in
  full_matrix.inst.in
  lapack_full_matrix.inst.in
  la_parallel_vector.inst.in
  la_parallel_block_vector.inst.in
  precondition_block.inst.in
  petsc_communication_pattern.inst.in
  relaxation_block.inst.in
  read_write_vector.inst.in
  solver.inst.in
  solver_gmres.inst.in
  sparse_matrix_ez.inst.in
  sparse_matrix.inst.in
  tensor_product_matrix.inst.in
  vector.inst.in
  vector_memory.inst.in
  vector_memory_release.inst.in
  )


# Add PETSc wrapper files. If PETSc has not been found,
# then these files should be empty and there is no need
# to even look at them
if(DEAL_II_WITH_PETSC)
  set(_unity_include_src
    ${_unity_include_src}
    petsc_communication_pattern.cc
    petsc_compatibility.cc
    petsc_full_matrix.cc
    petsc_matrix_base.cc
    petsc_matrix_free.cc
    petsc_parallel_block_sparse_matrix.cc
    petsc_parallel_block_vector.cc
    petsc_parallel_sparse_matrix.cc
    petsc_parallel_vector.cc
    petsc_precondition.cc
    petsc_snes.cc
    petsc_solver.cc
    petsc_ts.cc
    petsc_sparse_matrix.cc
    petsc_vector_base.cc
  )
endif()

# Same for SLEPc
if(DEAL_II_WITH_SLEPC)
  set(_unity_include_src
    ${_unity_include_src}
    slepc_solver.cc
    slepc_spectral_transformation.cc
  )
endif()

if(DEAL_II_WITH_GINKGO)
  set(_unity_include_src
    ${_unity_include_src}
    ginkgo_solver.cc
  )
endif()

# Also add Trilinos wrapper files
if(DEAL_II_WITH_TRILINOS)
  set(_unity_include_src
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
    trilinos_vector.cc
  )

  if(DEAL_II_TRILINOS_WITH_MUELU)
    set(_separate_src
      ${_separate_src}
      trilinos_precondition_muelu.cc
    )
  endif()

  if(DEAL_II_TRILINOS_WITH_TPETRA)
    set(_unity_include_src
      ${_unity_include_src}
      trilinos_tpetra_block_vector.cc
      trilinos_tpetra_block_sparse_matrix.cc
      trilinos_tpetra_communication_pattern.cc
      trilinos_tpetra_precondition.cc
      trilinos_tpetra_solver_direct.cc
      trilinos_tpetra_sparse_matrix.cc
      trilinos_tpetra_sparsity_pattern.cc
      trilinos_tpetra_vector.cc
    )
  endif()

  set(_inst
    ${_inst}
    trilinos_sparse_matrix.inst.in
    trilinos_vector.inst.in
  )
endif()

if(DEAL_II_WITH_SCALAPACK)
  set(_unity_include_src
    ${_unity_include_src}
    scalapack.cc
  )
  set(_inst
    ${_inst}
    scalapack.inst.in
  )
endif()


# determined by profiling
set(_n_includes_per_unity_file 20)

if(DEAL_II_UNITY_BUILD)
  # sort files so that the petsc and trilinos files are not together
  list(SORT _unity_include_src)
endif()

setup_source_list("${_unity_include_src}"
  "${_separate_src}"
  ${_n_includes_per_unity_file}
  _src
  )

file(GLOB _header CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/include/deal.II/lac/*.h
  )

define_object_library(object_lac OBJECT ${_src} ${_header} ${_inst})
expand_instantiations(object_lac "${_inst}")

collect_source_and_header_files("${_unity_include_src};${_separate_src}" "${_header}")
