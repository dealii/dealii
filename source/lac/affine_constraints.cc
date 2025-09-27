// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/affine_constraints.templates.h>


DEAL_II_NAMESPACE_OPEN

#include "lac/affine_constraints.inst"

/*
 * Note: You probably do not want to add your custom instantiation down
 * here but use affine_constraints.inst.in instead. We use the following
 * three macros for PETSc and Trilinos types because we lack a mechanism
 * for iterating over two "zipped" lists (we have to match trilinos/petsc
 * matrix and vector types).
 */

#define INSTANTIATE_DLTG_VECTOR(VectorType)                                    \
  template void AffineConstraints<VectorType::value_type>::condense<           \
    VectorType>(const VectorType &, VectorType &) const;                       \
  template void                                                                \
  AffineConstraints<VectorType::value_type>::condense<VectorType>(             \
    VectorType &) const;                                                       \
  template void                                                                \
  AffineConstraints<VectorType::value_type>::distribute_local_to_global<       \
    VectorType>(                                                               \
    const Vector<VectorType::value_type> &,                                    \
    const std::vector<AffineConstraints<VectorType::value_type>::size_type> &, \
    VectorType &,                                                              \
    const FullMatrix<VectorType::value_type> &) const;                         \
  template void                                                                \
  AffineConstraints<VectorType::value_type>::distribute_local_to_global<       \
    VectorType>(                                                               \
    const Vector<VectorType::value_type> &,                                    \
    const std::vector<AffineConstraints<VectorType::value_type>::size_type> &, \
    const std::vector<AffineConstraints<VectorType::value_type>::size_type> &, \
    VectorType &,                                                              \
    const FullMatrix<VectorType::value_type> &,                                \
    bool) const

#define INSTANTIATE_DLTG_VECTORMATRIX(MatrixType, VectorType) \
  template void AffineConstraints<MatrixType::value_type>::   \
    distribute_local_to_global<MatrixType, VectorType>(       \
      const FullMatrix<MatrixType::value_type> &,             \
      const Vector<VectorType::value_type> &,                 \
      const std::vector<AffineConstraints::size_type> &,      \
      MatrixType &,                                           \
      VectorType &,                                           \
      bool,                                                   \
      std::bool_constant<false>) const

#define INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(MatrixType, VectorType) \
  template void AffineConstraints<MatrixType::value_type>::         \
    distribute_local_to_global<MatrixType, VectorType>(             \
      const FullMatrix<MatrixType::value_type> &,                   \
      const Vector<VectorType::value_type> &,                       \
      const std::vector<AffineConstraints::size_type> &,            \
      MatrixType &,                                                 \
      VectorType &,                                                 \
      bool,                                                         \
      std::bool_constant<true>) const

#define INSTANTIATE_DLTG_MATRIX(MatrixType)                              \
  template void                                                          \
  AffineConstraints<MatrixType::value_type>::distribute_local_to_global< \
    MatrixType>(const FullMatrix<MatrixType::value_type> &,              \
                const std::vector<AffineConstraints::size_type> &,       \
                const std::vector<AffineConstraints::size_type> &,       \
                MatrixType &) const;                                     \
  template void                                                          \
  AffineConstraints<MatrixType::value_type>::distribute_local_to_global< \
    MatrixType>(const FullMatrix<MatrixType::value_type> &,              \
                const std::vector<AffineConstraints::size_type> &,       \
                const AffineConstraints<MatrixType::value_type> &,       \
                const std::vector<AffineConstraints::size_type> &,       \
                MatrixType &) const

#ifdef DEAL_II_WITH_PETSC
INSTANTIATE_DLTG_VECTOR(PETScWrappers::MPI::Vector);
INSTANTIATE_DLTG_VECTOR(PETScWrappers::MPI::BlockVector);

INSTANTIATE_DLTG_VECTORMATRIX(PETScWrappers::SparseMatrix, Vector<PetscScalar>);
INSTANTIATE_DLTG_VECTORMATRIX(PETScWrappers::SparseMatrix,
                              PETScWrappers::MPI::Vector);
INSTANTIATE_DLTG_VECTORMATRIX(PETScWrappers::MPI::SparseMatrix,
                              Vector<PetscScalar>);
INSTANTIATE_DLTG_VECTORMATRIX(PETScWrappers::MPI::SparseMatrix,
                              PETScWrappers::MPI::Vector);

INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(PETScWrappers::MPI::BlockSparseMatrix,
                                    Vector<PetscScalar>);
INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(PETScWrappers::MPI::BlockSparseMatrix,
                                    PETScWrappers::MPI::BlockVector);

INSTANTIATE_DLTG_MATRIX(PETScWrappers::SparseMatrix);
INSTANTIATE_DLTG_MATRIX(PETScWrappers::MPI::SparseMatrix);
INSTANTIATE_DLTG_MATRIX(PETScWrappers::MPI::BlockSparseMatrix);
#  ifndef DOXYGEN
#    ifdef DEAL_II_PETSC_WITH_COMPLEX
template void
dealii::AffineConstraints<double>::distribute<
  dealii::PETScWrappers::MPI::Vector>(
  dealii::PETScWrappers::MPI::Vector &) const;
#    endif
#  endif
#endif

#ifdef DEAL_II_WITH_TRILINOS
INSTANTIATE_DLTG_VECTOR(TrilinosWrappers::MPI::Vector);

INSTANTIATE_DLTG_VECTORMATRIX(TrilinosWrappers::SparseMatrix, Vector<double>);
INSTANTIATE_DLTG_VECTORMATRIX(TrilinosWrappers::SparseMatrix,
                              LinearAlgebra::distributed::Vector<double>);
INSTANTIATE_DLTG_VECTORMATRIX(TrilinosWrappers::SparseMatrix,
                              TrilinosWrappers::MPI::Vector);

INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(TrilinosWrappers::BlockSparseMatrix,
                                    Vector<double>);
INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(
  TrilinosWrappers::BlockSparseMatrix,
  LinearAlgebra::distributed::BlockVector<double>);
INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(TrilinosWrappers::BlockSparseMatrix,
                                    TrilinosWrappers::MPI::BlockVector);

INSTANTIATE_DLTG_MATRIX(TrilinosWrappers::SparseMatrix);
INSTANTIATE_DLTG_MATRIX(TrilinosWrappers::BlockSparseMatrix);

#  ifndef DOXYGEN
#    if defined(DEAL_II_TRILINOS_WITH_TPETRA) && defined(HAVE_TPETRA_INST_FLOAT)
// FIXME: This mixed variant is needed for multigrid and matrix free.
template void
dealii::AffineConstraints<double>::distribute<
  dealii::LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Host>>(
  dealii::LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Host> &)
  const;
template void
dealii::AffineConstraints<double>::distribute<
  dealii::LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Default>>(
  dealii::LinearAlgebra::TpetraWrappers::Vector<float, MemorySpace::Default> &)
  const;
#    endif
#  endif
#endif

#ifndef DOXYGEN
namespace internal
{
  namespace AffineConstraintsImplementation
  {
    template void
    set_zero_all(
      const std::vector<types::global_dof_index>                      &cm,
      LinearAlgebra::distributed::Vector<float, MemorySpace::Default> &vec);

    template void
    set_zero_all(
      const std::vector<types::global_dof_index>                       &cm,
      LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &vec);
  } // namespace AffineConstraintsImplementation
} // namespace internal

template void
AffineConstraints<float>::set_zero<
  LinearAlgebra::distributed::Vector<float, MemorySpace::Default>>(
  LinearAlgebra::distributed::Vector<float, MemorySpace::Default> &) const;

template void
AffineConstraints<double>::set_zero<
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default>>(
  LinearAlgebra::distributed::Vector<double, MemorySpace::Default> &) const;
#endif

DEAL_II_NAMESPACE_CLOSE
