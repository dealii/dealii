// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check the path of SolverCG that can apply loops to sub-ranges of the matrix
// in terms of the final solution vector, terminating at arbitrary numbers of
// iterations also when the solution has not converged yet, and make sure we
// obtain the same result as without interleaving.


#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>

#include "../tests.h"


struct MyDiagonalMatrix
{
  MyDiagonalMatrix(const LinearAlgebra::distributed::Vector<double> &diagonal)
    : diagonal(diagonal)
  {}

  void
  vmult(LinearAlgebra::distributed::Vector<double>       &dst,
        const LinearAlgebra::distributed::Vector<double> &src) const
  {
    dst = src;
    dst.scale(diagonal);
  }

  void
  vmult(
    LinearAlgebra::distributed::Vector<double>                        &dst,
    const LinearAlgebra::distributed::Vector<double>                  &src,
    const std::function<void(const unsigned int, const unsigned int)> &before,
    const std::function<void(const unsigned int, const unsigned int)> &after)
    const
  {
    before(0, dst.size());
    vmult(dst, src);
    after(0, dst.size());
  }

  const LinearAlgebra::distributed::Vector<double> &diagonal;
};



int
main()
{
  initlog();

  // Create diagonal matrix with entries between 1 and 15
  DiagonalMatrix<LinearAlgebra::distributed::Vector<double>> unit_matrix;
  unit_matrix.get_vector().reinit(15);
  unit_matrix.get_vector() = 1.0;

  LinearAlgebra::distributed::Vector<double> matrix_entries(unit_matrix.m());
  for (unsigned int i = 0; i < unit_matrix.m(); ++i)
    matrix_entries(i) = i + 1;
  MyDiagonalMatrix matrix(matrix_entries);

  LinearAlgebra::distributed::Vector<double> rhs(unit_matrix.m()),
    sol(unit_matrix.m());
  rhs = 1.;

  for (unsigned int n_iter = 1; n_iter < 12; ++n_iter)
    {
      deallog << "Test solution accuracy with " << n_iter << " iterations"
              << std::endl;
      deallog << "Solve with PreconditionIdentity: " << std::endl;
      IterationNumberControl control(n_iter, 1e-12);
      SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
      sol = 0;
      solver.solve(matrix, sol, rhs, PreconditionIdentity());
      sol.print(deallog.get_file_stream());

      deallog << "Solve with diagonal preconditioner: " << std::endl;
      sol = 0;
      solver.solve(matrix, sol, rhs, unit_matrix);
      sol.print(deallog.get_file_stream());
    }
}
