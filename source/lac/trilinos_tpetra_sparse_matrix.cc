// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/trilinos_tpetra_sparse_matrix.templates.h>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
// explicit instantiations
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
#    if (DEAL_II_TRILINOS_WITH_TPETRA_INST_DOUBLE)
    template class SparseMatrix<double, MemorySpace::Host>;

    template void
    SparseMatrix<double, MemorySpace::Host>::reinit(
      const IndexSet                       &parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<double, MemorySpace::Host>::reinit(
      const IndexSet                       &row_parallel_partitioning,
      const IndexSet                       &col_parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<double, MemorySpace::Host>::reinit(
      const dealii::DynamicSparsityPattern &);

    template void
    SparseMatrix<double, MemorySpace::Host>::vmult(
      Vector<double, MemorySpace::Host>       &dst,
      const Vector<double, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Host>::Tvmult(
      Vector<double, MemorySpace::Host>       &dst,
      const Vector<double, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Host>::vmult_add(
      Vector<double, MemorySpace::Host>       &dst,
      const Vector<double, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Host>::Tvmult_add(
      Vector<double, MemorySpace::Host>       &dst,
      const Vector<double, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Host>::vmult(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Host>::Tvmult(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Host>::vmult_add(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Host>::Tvmult_add(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;

    template class SparseMatrix<double, MemorySpace::Default>;

    template void
    SparseMatrix<double, MemorySpace::Default>::reinit(
      const IndexSet                       &parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<double, MemorySpace::Default>::reinit(
      const IndexSet                       &row_parallel_partitioning,
      const IndexSet                       &col_parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<double, MemorySpace::Default>::reinit(
      const dealii::DynamicSparsityPattern &);

    template void
    SparseMatrix<double, MemorySpace::Default>::vmult(
      Vector<double, MemorySpace::Default>       &dst,
      const Vector<double, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Default>::Tvmult(
      Vector<double, MemorySpace::Default>       &dst,
      const Vector<double, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Default>::vmult_add(
      Vector<double, MemorySpace::Default>       &dst,
      const Vector<double, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Default>::Tvmult_add(
      Vector<double, MemorySpace::Default>       &dst,
      const Vector<double, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Default>::vmult(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Default>::Tvmult(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Default>::vmult_add(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;

    template void
    SparseMatrix<double, MemorySpace::Default>::Tvmult_add(
      ::dealii::Vector<double>       &dst,
      const ::dealii::Vector<double> &src) const;
#    endif

#    if (DEAL_II_TRILINOS_WITH_TPETRA_INST_FLOAT)
    template class SparseMatrix<float, MemorySpace::Host>;

    template void
    SparseMatrix<float, MemorySpace::Host>::reinit(
      const IndexSet                       &parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<float, MemorySpace::Host>::reinit(
      const IndexSet                       &row_parallel_partitioning,
      const IndexSet                       &col_parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<float, MemorySpace::Host>::reinit(
      const dealii::DynamicSparsityPattern &);

    template void
    SparseMatrix<float, MemorySpace::Host>::vmult(
      Vector<float, MemorySpace::Host>       &dst,
      const Vector<float, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Host>::Tvmult(
      Vector<float, MemorySpace::Host>       &dst,
      const Vector<float, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Host>::vmult_add(
      Vector<float, MemorySpace::Host>       &dst,
      const Vector<float, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Host>::Tvmult_add(
      Vector<float, MemorySpace::Host>       &dst,
      const Vector<float, MemorySpace::Host> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Host>::vmult(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Host>::Tvmult(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Host>::vmult_add(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Host>::Tvmult_add(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;

    template class SparseMatrix<float, MemorySpace::Default>;

    template void
    SparseMatrix<float, MemorySpace::Default>::reinit(
      const IndexSet                       &parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<float, MemorySpace::Default>::reinit(
      const IndexSet                       &row_parallel_partitioning,
      const IndexSet                       &col_parallel_partitioning,
      const dealii::DynamicSparsityPattern &sparsity_pattern,
      const MPI_Comm                        communicator,
      const bool                            exchange_data);

    template void
    SparseMatrix<float, MemorySpace::Default>::reinit(
      const dealii::DynamicSparsityPattern &);

    template void
    SparseMatrix<float, MemorySpace::Default>::vmult(
      Vector<float, MemorySpace::Default>       &dst,
      const Vector<float, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Default>::Tvmult(
      Vector<float, MemorySpace::Default>       &dst,
      const Vector<float, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Default>::vmult_add(
      Vector<float, MemorySpace::Default>       &dst,
      const Vector<float, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Default>::Tvmult_add(
      Vector<float, MemorySpace::Default>       &dst,
      const Vector<float, MemorySpace::Default> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Default>::vmult(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Default>::Tvmult(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Default>::vmult_add(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;

    template void
    SparseMatrix<float, MemorySpace::Default>::Tvmult_add(
      ::dealii::Vector<float>       &dst,
      const ::dealii::Vector<float> &src) const;
#    endif
  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
