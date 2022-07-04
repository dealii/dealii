// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2022 by the deal.II authors
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

#ifndef dealii_mg_transfer_matrix_free_templates_h
#define dealii_mg_transfer_matrix_free_templates_h

#include <deal.II/base/config.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>


DEAL_II_NAMESPACE_OPEN

template <int dim, typename Number, typename TransferType>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::prolongate(
  const unsigned int                                     to_level,
  LinearAlgebra::distributed::BlockVector<Number> &      dst,
  const LinearAlgebra::distributed::BlockVector<Number> &src) const
{
  const unsigned int n_blocks = src.n_blocks();
  AssertDimension(dst.n_blocks(), n_blocks);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      const unsigned int data_block = same_for_all ? 0 : b;
      get_matrix_free_transfer(data_block)
        .prolongate(to_level, dst.block(b), src.block(b));
    }
}



template <int dim, typename Number, typename TransferType>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::prolongate_and_add(
  const unsigned int                                     to_level,
  LinearAlgebra::distributed::BlockVector<Number> &      dst,
  const LinearAlgebra::distributed::BlockVector<Number> &src) const
{
  const unsigned int n_blocks = src.n_blocks();
  AssertDimension(dst.n_blocks(), n_blocks);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      const unsigned int data_block = same_for_all ? 0 : b;
      get_matrix_free_transfer(data_block)
        .prolongate_and_add(to_level, dst.block(b), src.block(b));
    }
}



template <int dim, typename Number, typename TransferType>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::restrict_and_add(
  const unsigned int                                     from_level,
  LinearAlgebra::distributed::BlockVector<Number> &      dst,
  const LinearAlgebra::distributed::BlockVector<Number> &src) const
{
  const unsigned int n_blocks = src.n_blocks();
  AssertDimension(dst.n_blocks(), n_blocks);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      const unsigned int data_block = same_for_all ? 0 : b;
      get_matrix_free_transfer(data_block)
        .restrict_and_add(from_level, dst.block(b), src.block(b));
    }
}

DEAL_II_NAMESPACE_CLOSE

#endif
