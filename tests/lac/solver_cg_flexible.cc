// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check flexible variant of SolverCG by verifying that it converges more
// quickly in case of a variable preconditioner


#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


// Class to represent an approximate inverse with variable number of
// iterations to stress the flexible CG implementation
class IterativeInverse
{
public:
  IterativeInverse(const DiagonalMatrix<Vector<double>> &matrix)
    : matrix(matrix)
    , apply_count(0)
  {}

  void
  vmult(Vector<double> &dst, const Vector<double> &src) const
  {
    IterationNumberControl control(apply_count % 2 ? 4 : 3, 1e-12);
    SolverCG<>             inner_solver(control);

    dst = 0;
    inner_solver.solve(matrix, dst, src, PreconditionIdentity());

    ++apply_count;
  }

private:
  const DiagonalMatrix<Vector<double>> &matrix;
  mutable unsigned int                  apply_count;
};



SolverControl::State
monitor_norm(const unsigned int iteration,
             const double       check_value,
             const Vector<double> &)
{
  deallog << "   CG estimated residual at iteration " << iteration << ": "
          << check_value << std::endl;
  return SolverControl::success;
}


int
main()
{
  initlog();

  // Create diagonal matrix with entries between 1 and 30
  DiagonalMatrix<Vector<double>> matrix;
  matrix.get_vector().reinit(30);
  for (unsigned int i = 0; i < matrix.m(); ++i)
    matrix.get_vector()[i] = i + 1.0;

  Vector<double> rhs(matrix.m()), sol(matrix.m());
  rhs = 1.;

  IterativeInverse variable_preconditioner(matrix);

  // Try to solve with standard CG method - note that we solve with a coarse
  // tolerance to avoid trouble with roundoff
  {
    SolverControl control(30, 1e-4);
    SolverCG<>    solver(control);
    solver.connect(&monitor_norm);
    solver.solve(matrix, sol, rhs, variable_preconditioner);
  }

  // And now use the flexible CG method
  {
    SolverControl      control(30, 1e-4);
    SolverFlexibleCG<> solver(control);
    solver.connect(&monitor_norm);
    sol = 0;
    solver.solve(matrix, sol, rhs, variable_preconditioner);
  }
}
