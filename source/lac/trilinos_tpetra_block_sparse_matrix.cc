// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/config.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA

#  include <deal.II/lac/trilinos_tpetra_block_sparse_matrix.templates.h>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
// explicit instantiations
namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template class BlockSparseMatrix<double, MemorySpace::Host>;

    template void
    BlockSparseMatrix<double, MemorySpace::Host>::reinit(
      const ::dealii::BlockDynamicSparsityPattern &);

    template void
    BlockSparseMatrix<double, MemorySpace::Host>::reinit<
      ::dealii::BlockDynamicSparsityPattern>(
      const std::vector<dealii::IndexSet> &,
      const ::dealii::BlockDynamicSparsityPattern &,
      MPI_Comm,
      bool);

    template void
    BlockSparseMatrix<double, MemorySpace::Host>::vmult(
      TpetraWrappers::Vector<double, MemorySpace::Host> &,
      const TpetraWrappers::Vector<double, MemorySpace::Host> &) const;

    template void
    BlockSparseMatrix<double, MemorySpace::Host>::Tvmult(
      TpetraWrappers::Vector<double, MemorySpace::Host> &,
      const TpetraWrappers::Vector<double, MemorySpace::Host> &) const;

    template void
    BlockSparseMatrix<double, MemorySpace::Host>::vmult(
      TpetraWrappers::BlockVector<double, MemorySpace::Host> &,
      const TpetraWrappers::BlockVector<double, MemorySpace::Host> &) const;

    template void
    BlockSparseMatrix<double, MemorySpace::Host>::Tvmult(
      TpetraWrappers::BlockVector<double, MemorySpace::Host> &,
      const TpetraWrappers::BlockVector<double, MemorySpace::Host> &) const;
  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
