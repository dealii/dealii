// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/affine_constraints.templates.h>


DEAL_II_NAMESPACE_OPEN

#include "affine_constraints.inst"

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
      std::integral_constant<bool, false>) const

#define INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(MatrixType, VectorType) \
  template void AffineConstraints<MatrixType::value_type>::         \
    distribute_local_to_global<MatrixType, VectorType>(             \
      const FullMatrix<MatrixType::value_type> &,                   \
      const Vector<VectorType::value_type> &,                       \
      const std::vector<AffineConstraints::size_type> &,            \
      MatrixType &,                                                 \
      VectorType &,                                                 \
      bool,                                                         \
      std::integral_constant<bool, true>) const

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
#endif

#ifdef DEAL_II_WITH_TRILINOS
INSTANTIATE_DLTG_VECTOR(TrilinosWrappers::MPI::Vector);

INSTANTIATE_DLTG_VECTORMATRIX(TrilinosWrappers::SparseMatrix, Vector<double>);
INSTANTIATE_DLTG_VECTORMATRIX(TrilinosWrappers::SparseMatrix,
                              TrilinosWrappers::MPI::Vector);

INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(TrilinosWrappers::BlockSparseMatrix,
                                    Vector<double>);
INSTANTIATE_DLTG_BLOCK_VECTORMATRIX(TrilinosWrappers::BlockSparseMatrix,
                                    TrilinosWrappers::MPI::BlockVector);

INSTANTIATE_DLTG_MATRIX(TrilinosWrappers::SparseMatrix);
INSTANTIATE_DLTG_MATRIX(TrilinosWrappers::BlockSparseMatrix);
#endif

/*
 * Allocate scratch data.
 *
 * We cannot use the generic template instantiation because we need to
 * provide an initializer object of type
 * internals::AffineConstraintsData<Number> that can be passed to the
 * constructor of scratch_data (it won't allow one to be constructed in
 * place).
 */

namespace internal
{
  namespace AffineConstraints
  {
#define SCRATCH_INITIALIZER(number, Name)              \
  ScratchData<number> scratch_data_initializer_##Name; \
  template <>                                          \
  Threads::ThreadLocalStorage<ScratchData<number>>     \
    AffineConstraintsData<number>::scratch_data(       \
      scratch_data_initializer_##Name)

    SCRATCH_INITIALIZER(double, d);
    SCRATCH_INITIALIZER(float, f);
#ifdef DEAL_II_WITH_COMPLEX_VALUES
    SCRATCH_INITIALIZER(std::complex<double>, cd);
    SCRATCH_INITIALIZER(std::complex<float>, cf);
#endif
#undef SCRATCH_INITIALIZER
  } // namespace AffineConstraints
} // namespace internal

DEAL_II_NAMESPACE_CLOSE
