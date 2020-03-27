// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Test non symmetric variants:

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#define PRINTME(name, var)                              \
  deallog << "Block vector: " name << ":" << std::endl; \
  for (unsigned int i = 0; i < var.n_blocks(); ++i)     \
    deallog << "[block " << i << " ]  " << var.block(i);



int
main()
{
  initlog();
  deallog << std::setprecision(10);

  // SparseMatrix:
  {
    SparsityPattern sparsity_pattern(10, 5, 0);
    sparsity_pattern.compress();

    SparseMatrix<double> a(sparsity_pattern);

    auto op_a = linear_operator(a);

    Vector<double> u;
    op_a.reinit_domain_vector(u, false);
    Vector<double> v;
    op_a.reinit_range_vector(v, false);
    op_a.vmult(v, u);
    op_a.vmult_add(v, u);
    op_a.Tvmult(u, v);
    op_a.Tvmult_add(u, v);

    deallog << "OK" << std::endl;
  }

  // BlockSparseMatrix:
  {
    BlockDynamicSparsityPattern dsp(5, 3);
    for (unsigned int i = 0; i < 5; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        dsp.block(i, j).reinit(10, 5);
    dsp.collect_sizes();

    BlockSparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.compress();

    BlockSparseMatrix<double> a(sparsity_pattern);

    auto op_a = linear_operator<BlockVector<double>>(a);

    BlockVector<double> u;
    op_a.reinit_domain_vector(u, false);
    BlockVector<double> v;
    op_a.reinit_range_vector(v, false);
    op_a.vmult(v, u);
    op_a.vmult_add(v, u);
    op_a.Tvmult(u, v);
    op_a.Tvmult_add(u, v);
    deallog << "OK" << std::endl;
  }
}
