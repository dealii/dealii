// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test block_operator's constructors and copy assignment variants

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


int
main()
{
  initlog();

  {
    BlockSparseMatrix<double> a;

    BlockLinearOperator<> op_b(a);
    auto                  op_c = block_operator<BlockVector<double>>(a);
    op_c                       = a;

    auto op_d = block_operator(a);
    op_d      = a;
  }

  auto op_a = LinearOperator<>();


  {
    std::array<decltype(op_a), 2>                a{{op_a, op_a}};
    std::array<std::array<decltype(op_a), 2>, 2> temp{{a, a}};

    BlockLinearOperator<> op_b(temp);
    auto op_c = block_operator<2, 2, BlockVector<double>>(temp);
    op_c      = temp;

    auto op_d = block_operator(temp);
    op_d      = temp;
  }

  {
    std::array<decltype(op_a), 2> temp{{op_a, op_a}};

    auto op_c = block_diagonal_operator<2, BlockVector<double>>(temp);
    op_c      = temp;

    auto op_d = block_diagonal_operator(temp);
    op_d      = temp;
  }

  deallog << "OK!" << std::endl;
}
