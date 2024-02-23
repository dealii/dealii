// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2024 by the deal.II authors
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
    template class BlockSparseMatrix<double>;

    template void
    BlockSparseMatrix<double>::reinit(
      const ::dealii::BlockDynamicSparsityPattern &);

    template void
    BlockSparseMatrix<double>::vmult(
      TpetraWrappers::Vector<double> &,
      const TpetraWrappers::Vector<double> &) const;

    template void
    BlockSparseMatrix<double>::Tvmult(
      TpetraWrappers::Vector<double> &,
      const TpetraWrappers::Vector<double> &) const;

    template void
    BlockSparseMatrix<double>::vmult(
      TpetraWrappers::BlockVector<double> &,
      const TpetraWrappers::BlockVector<double> &) const;

    template void
    BlockSparseMatrix<double>::Tvmult(
      TpetraWrappers::BlockVector<double> &,
      const TpetraWrappers::BlockVector<double> &) const;
  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
