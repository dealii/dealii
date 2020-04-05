// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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


// adapted from a testcase by Roger Young, sent to the mailing list
// 2007-03-02, that illustrates that bicgstab can't handle early
// success

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  GrowingVectorMemory<> mem;
  SolverControl         control(100, 1.e-3);
  SolverBicgstab<>      bicgstab(control, mem);

  SparsityPattern sparsity_pattern(4, 4, 4);
  for (unsigned int i = 0; i < 4; ++i)
    for (unsigned int j = 0; j < 4; ++j)
      sparsity_pattern.add(i, j);
  sparsity_pattern.compress();

  SparseMatrix<double> M(sparsity_pattern);
  M.set(0, 0, 21.1);
  M.set(0, 1, 0);
  M.set(0, 2, 0);
  M.set(0, 3, 0);
  M.set(1, 1, 7.033333333);
  M.set(1, 0, 0);
  M.set(1, 2, 0);
  M.set(1, 3, 3.516666667);
  M.set(2, 2, 21.1);
  M.set(2, 0, 0);
  M.set(2, 1, 0);
  M.set(2, 3, 0);
  M.set(3, 3, 7.033333333);
  M.set(3, 0, 0);
  M.set(3, 1, 3.516666667);
  M.set(3, 2, 0);

  Vector<double> rhs(4);
  rhs(0) = rhs(2) = 0;
  rhs(1) = rhs(3) = 0.0975;


  SparseILU<double>::AdditionalData data;
  data.use_this_sparsity = &sparsity_pattern;
  SparseILU<double> ilu;
  ilu.initialize(M, data);

  Vector<double> solution(4);
  // this would fail with elements of
  // the solution vector being set to
  // NaN before Roger's suggested
  // change
  bicgstab.solve(M, solution, rhs, ilu);

  for (unsigned int i = 0; i < 4; ++i)
    deallog << solution(i) << std::endl;

  Vector<double> res(4);
  M.residual(res, solution, rhs);
  deallog << "residual=" << res.l2_norm() << std::endl;
}
