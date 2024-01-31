// ---------------------------------------------------------------------
//
// Copyright (C) 2024 by the deal.II authors
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

    template void
    BlockSparseMatrix<double>::vmult(
      ::dealii::BlockVector<double> &,
      const ::dealii::BlockVector<double> &) const;

    template void
    BlockSparseMatrix<double>::Tvmult(
      ::dealii::BlockVector<double> &,
      const ::dealii::BlockVector<double> &) const;


  } // namespace TpetraWrappers
} // namespace LinearAlgebra
#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_TRILINOS_WITH_TPETRA
