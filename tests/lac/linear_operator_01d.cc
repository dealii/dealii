// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2019 by the deal.II authors
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

// Verify that the inverse_operator function signature without a
// preconditioner argument and the one with PreconditionIdentity() as a
// temporary works correctly.

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


int
main()
{
  initlog();

  Vector<double> answer(2);
  answer[0] = 0.25;
  answer[1] = 0.25;

  SparsityPattern sparsity_pattern(2, 2, 1);
  sparsity_pattern.add(0, 0);
  sparsity_pattern.add(1, 1);
  sparsity_pattern.compress();

  SparseMatrix<double> A(sparsity_pattern);
  A.set(0, 0, 4.0);
  A.set(1, 1, 4.0);

  Vector<double> b(2);
  b[0] = 1.0;
  b[1] = 1.0;

  const auto lo_A = linear_operator(A);

  // As we all remember from numerical analysis class, CG will converge in
  // at most two iterations
  SolverControl solver_control_A(2, 1.0e-15);

  SolverCG<Vector<double>> solver_A(solver_control_A);

  const auto lo_A_inv   = inverse_operator(lo_A, solver_A);
  const auto lo_A_inv_t = transpose_operator(lo_A_inv);

  deallog.depth_file(0);

  unsigned int n_mistakes{0};

  Vector<double> residual; // keep storage location to trigger bug.

  for (unsigned int j = 0; j < 1000; ++j)
    {
      // test Tvmult:
      residual = lo_A_inv_t * b;
      residual -= answer;
      if (residual.l2_norm() > 1e-10)
        ++n_mistakes;

      // test Tvmult_add:
      residual = lo_A_inv_t * b - answer;
      if (residual.l2_norm() > 1e-10)
        ++n_mistakes;
    }

  deallog.depth_file(3);

  deallog << "number of mistakes: " << n_mistakes << std::endl;
  if (n_mistakes != 0)
    deallog << "Ooops!" << std::endl;
}
