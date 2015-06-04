// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Test upper_triangular_block_op and lower_triangular_operator
// simultaneously:

#include "../tests.h"

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#define PRINTME(name, var) \
  deallog << "Block vector: " name << ":" << std::endl; \
  for (unsigned int i = 0; i < var.n_blocks(); ++i) \
    deallog << "[block " << i << " ]  " << var.block(i);


using namespace dealii;

int main()
{
  initlog();
  deallog << std::setprecision(2);

  // BlockSparseMatrix:
  {
    BlockDynamicSparsityPattern dsp(3, 3);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        dsp.block(i, j).reinit (1, 1);
    dsp.collect_sizes ();

    BlockSparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from(dsp);
    sparsity_pattern.compress();

    BlockSparseMatrix<double> a (sparsity_pattern);
    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
        a.block(i,j).set(0, 0, 2 + j + i);

    auto op_a = linear_operator<BlockVector<double>>(a);

    auto lower_triangular_block_op = lower_triangular_operator<3, BlockVector<double>, BlockVector<double>, BlockSparseMatrix<double> >(a);
    auto upper_triangular_block_op = upper_triangular_operator<3, BlockVector<double>, BlockVector<double>, BlockSparseMatrix<double> >(a);


    BlockVector<double> u1;
    op_a.reinit_domain_vector(u1, false);
    BlockVector<double> v1;
    op_a.reinit_range_vector(v1, false);
    BlockVector<double> u2;
    op_a.reinit_domain_vector(u2, false);
    BlockVector<double> v2;
    op_a.reinit_range_vector(v2, false);


    // check 1
    for(unsigned int i = 0; i<3; ++i)
    {
      u1.block(i)[0] = i+1;
      u2.block(i)[0] = i+1;
      v1.block(i)[0] = 1;
      v2.block(i)[0] = 1;
    }

    lower_triangular_block_op.vmult(v1, u1);
    upper_triangular_block_op.Tvmult(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);


    // check 2
    for(unsigned int i = 0; i<3; ++i)
    {
      u1.block(i)[0] = i+1;
      u2.block(i)[0] = i+1;
      v1.block(i)[0] = 1;
      v2.block(i)[0] = 1;
    }

    lower_triangular_block_op.Tvmult(v1, u1);
    upper_triangular_block_op.vmult(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);

    // check 3
    for(unsigned int i = 0; i<3; ++i)
    {
      u1.block(i)[0] = i+1;
      u2.block(i)[0] = i+1;
      v1.block(i)[0] = 1;
      v2.block(i)[0] = 1;
    }

    lower_triangular_block_op.vmult_add(v1, u1);
    upper_triangular_block_op.Tvmult_add(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);


    // check 4
    for(unsigned int i = 0; i<3; ++i)
    {
      u1.block(i)[0] = i+1;
      u2.block(i)[0] = i+1;
      v1.block(i)[0] = 1;
      v2.block(i)[0] = 1;
    }

    lower_triangular_block_op.Tvmult_add(v1, u1);
    upper_triangular_block_op.vmult_add(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);
  }
}
