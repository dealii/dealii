// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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


// test the umfpack sparse direct solver on a simple 2x2 block matrix
// that equals the unit matrix. same as umfpack_05, but with the
// difference that there, the matrix is partitioned into 2+2 rows and
// columns. here, do it as 3+1

#include "../tests.h"
#include <iostream>
#include <fstream>
#include <cstdlib>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_direct.h>


void test ()
{
  std::vector<unsigned int> size (2);
  size[0] = 3;
  size[1] = 1;

  BlockSparsityPattern b_sparsity_pattern;
  b_sparsity_pattern.reinit(size.size(),size.size());
  for (unsigned int k=0; k<size.size(); ++k)
    for (unsigned int l=0; l<size.size(); ++l)
      b_sparsity_pattern.block(k,l).reinit(size[k],
                                           size[l],
                                           3);
  b_sparsity_pattern.collect_sizes();
  for (unsigned int i=0; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      b_sparsity_pattern.add (i,j);
  b_sparsity_pattern.compress();

  BlockSparseMatrix<double> Bb (b_sparsity_pattern);
  for (unsigned int i=0; i<4; ++i)
    Bb.add (i,i,1);

  SparseDirectUMFPACK umfpackb;
  umfpackb.factorize(Bb);

  Vector<double> ubb(4);
  for (unsigned int i=0; i<4; ++i)
    ubb(i) = i;

  umfpackb.solve(ubb);

  for (unsigned int i=0; i<4; ++i)
    Assert (std::fabs (ubb(i) - i) < 1e-12,
            ExcInternalError());

  deallog << "OK" << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();
}
