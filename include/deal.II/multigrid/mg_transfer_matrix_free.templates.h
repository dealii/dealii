// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_transfer_matrix_free_templates_h
#define dealii_mg_transfer_matrix_free_templates_h

#include <deal.II/base/config.h>

#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <type_traits>


DEAL_II_NAMESPACE_OPEN

template <int dim, typename Number, typename TransferType>
void
MGTransferBlockMatrixFreeBase<dim, Number, TransferType>::prolongate(
  const unsigned int                                     to_level,
  LinearAlgebra::distributed::BlockVector<Number>       &dst,
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
  LinearAlgebra::distributed::BlockVector<Number>       &dst,
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
  LinearAlgebra::distributed::BlockVector<Number>       &dst,
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
