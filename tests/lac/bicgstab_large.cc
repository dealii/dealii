// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2018 by the deal.II authors
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


// check that bicgstab does not exit early when very large matrices are used

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"

#include "../testmatrix.h"


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  SparsityPattern sparsity_pattern(4, 4);
  sparsity_pattern.compress();

  SparseMatrix<double> M(sparsity_pattern);
  M.diag_element(0) = 1;
  M.diag_element(1) = 10;
  M.diag_element(2) = 11;
  M.diag_element(3) = 42;

  Vector<double> rhs(4);
  rhs = 1;

  Vector<double> solution(4);

  {
    SolverControl    control(100, 1.e-3);
    SolverBicgstab<> bicgstab(control);
    bicgstab.solve(M, solution, rhs, PreconditionIdentity());
  }

  solution.print(deallog.get_file_stream());

  Vector<double> res(4);
  M.residual(res, solution, rhs);
  deallog << "residual=" << res.l2_norm() << std::endl;

  // now set up the same problem but with matrix entries scaled by 1e10 and
  // solver tolerance scaled by 1e10. should get the same solution
  SparseMatrix<double> M1(sparsity_pattern);
  M1.add(1e10, M);
  rhs *= 1e10;
  solution = 0;

  {
    SolverControl    control(100, 1.e7);
    SolverBicgstab<> bicgstab(control);
    bicgstab.solve(M1, solution, rhs, PreconditionIdentity());
  }
  solution.print(deallog.get_file_stream());
  M1.residual(res, solution, rhs);
  deallog << "residual=" << res.l2_norm() << std::endl;
}
