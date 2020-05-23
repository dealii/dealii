// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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


// MinRes can't deal with block systems at the time of writing this test


#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_minres.h>

#include "../tests.h"


int
main()
{
  initlog();
  deallog << std::setprecision(2);


  // assemble a 2x2 block identity
  // matrix
  BlockSparsityPattern block_structure(2, 2);
  block_structure.block(0, 0).reinit(2, 2, 1);
  block_structure.block(0, 1).reinit(2, 2, 1);
  block_structure.block(1, 0).reinit(2, 2, 1);
  block_structure.block(1, 1).reinit(2, 2, 1);
  block_structure.collect_sizes();

  block_structure.add(0, 0);
  block_structure.add(1, 1);
  block_structure.add(2, 2);
  block_structure.add(3, 3);
  block_structure.compress();

  BlockSparseMatrix<double> block_A(block_structure);
  block_A.add(0, 0, 1);
  block_A.add(1, 1, 1);
  block_A.add(2, 2, 1);
  block_A.add(3, 3, 1);

  // block vector with 2 blocks of 2
  // components each
  BlockVector<double> a(2, 2);
  a(0) = 2;
  a(1) = 3;
  a(2) = 4;
  a(3) = 5;

  // this should work (check that
  // sizes of objects are all
  // correct)
  BlockVector<double> b(2, 2);
  block_A.vmult(b, a);

  // and while at it also check
  // whether the result was correct
  AssertThrow(b(0) == 2, ExcInternalError());
  AssertThrow(b(1) == 3, ExcInternalError());
  AssertThrow(b(2) == 4, ExcInternalError());
  AssertThrow(b(3) == 5, ExcInternalError());

  // now solve with MinRes. This
  // didn't work at one point
  SolverControl                     solver_control(1000, 1e-12);
  SolverMinRes<BlockVector<double>> minres(solver_control);

  minres.solve(block_A, b, a, PreconditionIdentity());

  deallog << "OK" << std::endl;
}
