// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test block_back_substitution and block_forward_substitution:

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#define PRINTME(name, var)                                            \
  deallog << "Block vector: " name << ':' << std::endl;               \
  for (unsigned int i = 0; i < var.n_blocks(); ++i)                   \
    deallog << "[block " << i << " ]  " << var.block(i) << std::endl; \
  deallog << std::endl;



int
main()
{
  initlog();
  deallog << std::setprecision(12);

  // BlockSparseMatrix:
  {
    BlockDynamicSparsityPattern dsp(3, 3);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        dsp.block(i, j).reinit(1, 1);
    dsp.collect_sizes();

    BlockSparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.compress();

    BlockSparseMatrix<double> a(sparsity_pattern);
    for (unsigned int i = 0; i < 3; ++i)
      {
        a.block(i, i).set(0, 0, i + i + 1);
        for (unsigned int j = 0; j < i; ++j)
          a.block(i, j).set(0, 0, 10);
      }

    BlockSparseMatrix<double> d(sparsity_pattern);
    for (unsigned int i = 0; i < 3; ++i)
      d.block(i, i).set(0, 0, 1.0 / (i + i + 1));

    auto op_a = linear_operator<BlockVector<double>>(a);

    auto op_b1 = block_operator(a);
    auto op_b2 = block_diagonal_operator(d);

    auto inverse_op_a =
      block_forward_substitution<BlockVector<double>>(op_b1, op_b2);

    auto identity = inverse_op_a * op_a;

    BlockVector<double> u;
    BlockVector<double> v;

    deallog << " -- Matrix -- " << std::endl;
    op_a.reinit_domain_vector(u, false);
    op_a.reinit_range_vector(v, false);
    for (unsigned int j = 0; j < 3; ++j)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            u.block(i)[0] = 0;
            v.block(i)[0] = 0;
          }
        u.block(j)[0] = 1;

        op_a.vmult(v, u);

        PRINTME("v", v);
      }

    deallog << " -- Inverse -- " << std::endl;
    inverse_op_a.reinit_domain_vector(u, false);
    inverse_op_a.reinit_range_vector(v, true);
    for (unsigned int j = 0; j < 3; ++j)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            u.block(i)[0] = 0;
            v.block(i)[0] = 0;
          }
        u.block(j)[0] = 1;

        inverse_op_a.vmult(v, u);

        PRINTME("v", v);
      }

    deallog << " -- Identity -- " << std::endl;
    identity.reinit_domain_vector(u, false);
    identity.reinit_range_vector(v, false);
    for (unsigned int j = 0; j < 3; ++j)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            u.block(i)[0] = 0;
            v.block(i)[0] = 0;
          }
        u.block(j)[0] = 1;

        identity.vmult(v, u);

        PRINTME("v", v);
      }
  }


  {
    BlockDynamicSparsityPattern dsp(3, 3);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        dsp.block(i, j).reinit(1, 1);
    dsp.collect_sizes();

    BlockSparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.compress();

    BlockSparseMatrix<double> a(sparsity_pattern);
    for (unsigned int i = 0; i < 3; ++i)
      {
        a.block(i, i).set(0, 0, i + i + 1);
        for (unsigned int j = i + 1; j < 3; ++j)
          a.block(i, j).set(0, 0, 10);
      }

    BlockSparseMatrix<double> d(sparsity_pattern);
    for (unsigned int i = 0; i < 3; ++i)
      d.block(i, i).set(0, 0, 1.0 / (i + i + 1));

    auto op_a = linear_operator<BlockVector<double>>(a);

    auto op_b1 = block_operator(a);
    auto op_b2 = block_diagonal_operator(d);

    auto inverse_op_a =
      block_back_substitution<BlockVector<double>>(op_b1, op_b2);

    auto identity = inverse_op_a * op_a;

    BlockVector<double> u;
    BlockVector<double> v;

    deallog << " -- Matrix -- " << std::endl;
    op_a.reinit_domain_vector(u, false);
    op_a.reinit_range_vector(v, false);
    for (unsigned int j = 0; j < 3; ++j)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            u.block(i)[0] = 0;
            v.block(i)[0] = 0;
          }
        u.block(j)[0] = 1;

        op_a.vmult(v, u);

        PRINTME("v", v);
      }

    deallog << " -- Inverse -- " << std::endl;
    inverse_op_a.reinit_domain_vector(u, false);
    inverse_op_a.reinit_range_vector(v, true);
    for (unsigned int j = 0; j < 3; ++j)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            u.block(i)[0] = 0;
            v.block(i)[0] = 0;
          }
        u.block(j)[0] = 1;

        inverse_op_a.vmult(v, u);

        PRINTME("v", v);
      }

    deallog << " -- Identity -- " << std::endl;
    identity.reinit_domain_vector(u, false);
    identity.reinit_range_vector(v, false);
    for (unsigned int j = 0; j < 3; ++j)
      {
        for (unsigned int i = 0; i < 3; ++i)
          {
            u.block(i)[0] = 0;
            v.block(i)[0] = 0;
          }
        u.block(j)[0] = 1;

        identity.vmult(v, u);

        PRINTME("v", v);
      }
  }
}
