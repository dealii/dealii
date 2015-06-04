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
// (in the case of array of array of LinearOperator as argument)
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
  for (unsigned int j = 0; j < 2; ++j) \
      deallog << "[block " << j  << " ]  " << var.block(j);

void reinit_vec(BlockVector<double> &u1, BlockVector<double> &u2, BlockVector<double> &v1, BlockVector<double> &v2)
{
  for(unsigned int i = 0; i<2; ++i)
  {
    for(unsigned int j = 0; j<2; ++j)
    {
    u1.block(i)[j] = i+j;
    u2.block(i)[j] = i+j;
    v1.block(i)[j] = 1;
    v2.block(i)[j] = 1;
  }
  }
}

using namespace dealii;

int main()
{
  initlog();
  deallog << std::setprecision(2);

  // BlockSparseMatrix:
  {
    BlockSparsityPattern dsp(2,2);
    // set sizes
    for (unsigned int i=0; i<2; ++i)
      for (unsigned int j=0; j<2; ++j)
        dsp.block(i,j).reinit ( 2, 2, 5);
    dsp.collect_sizes ();

    for (unsigned int row=0; row<4; ++row)
      for (unsigned int i=0; i<4; ++i)
        dsp.add (row, i);
    dsp.compress ();

    BlockSparseMatrix<double> a (dsp);

    for (unsigned int i = 0; i < 2; ++i)
      for (unsigned int j = 0; j < 2; ++j)
        for (unsigned int h = 0; h < 2; ++h)
          for (unsigned int k = 0; k < 2; ++k)
            a.block(i,j).set(h, k, j + i + h + k);

    auto op00 = linear_operator<Vector<double>,Vector<double>>(a.block(0,0));
    auto op01 = linear_operator<Vector<double>,Vector<double>>(a.block(0,1));
    auto op10 = linear_operator<Vector<double>,Vector<double>>(a.block(1,0));
    auto op11 = linear_operator<Vector<double>,Vector<double>>(a.block(1,1));

    auto lower_triangular_block_op =
          lower_triangular_operator<  2,
                                      BlockVector<double>,
                                      BlockVector<double>>({{
      {{ op00, op01 }} ,
      {{ op10, op11 }}
    }
  });
    auto upper_triangular_block_op =
          upper_triangular_operator<  2,
                                      BlockVector<double>,
                                      BlockVector<double>>({{
      {{ op00, op01 }} ,
      {{ op10, op11 }}
    }
  });

    BlockVector<double> u1;
    lower_triangular_block_op.reinit_domain_vector(u1, false);
    BlockVector<double> v1;
    lower_triangular_block_op.reinit_range_vector(v1, false);
    BlockVector<double> u2;
    upper_triangular_block_op.reinit_domain_vector(u2, false);
    BlockVector<double> v2;
    upper_triangular_block_op.reinit_range_vector(v2, false);


    // check 1
    reinit_vec(u1,u2,v1,v2);
    PRINTME("u1", u1);
    PRINTME("u2", u2);
    PRINTME("v1", v1);
    PRINTME("v2", v2);

    lower_triangular_block_op.vmult(v1, u1);
    upper_triangular_block_op.Tvmult(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);


    // check 2
    reinit_vec(u1,u2,v1,v2);

    lower_triangular_block_op.Tvmult(v1, u1);
    upper_triangular_block_op.vmult(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);

    // check 3
    reinit_vec(u1,u2,v1,v2);

    lower_triangular_block_op.vmult_add(v1, u1);
    upper_triangular_block_op.Tvmult_add(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);


    // check 4
    reinit_vec(u1,u2,v1,v2);

    lower_triangular_block_op.Tvmult_add(v1, u1);
    upper_triangular_block_op.vmult_add(v2, u2);

    PRINTME("v1", v1);
    PRINTME("v2", v2);
  }
}
