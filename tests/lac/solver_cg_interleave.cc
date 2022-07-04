// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2022 by the deal.II authors
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

// Check the path of SolverCG that can apply loops to sub-ranges of the
// matrix. That feature is used for matrix-free loops with increased data
// locality.


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
  vmult(LinearAlgebra::distributed::Vector<double> &      dst,
        const LinearAlgebra::distributed::Vector<double> &src) const
  {
    dst = src;
    dst.scale(diagonal);
  }

  void
  vmult(
    LinearAlgebra::distributed::Vector<double> &                       dst,
    const LinearAlgebra::distributed::Vector<double> &                 src,
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



SolverControl::State
monitor_norm(const unsigned int iteration,
             const double       check_value,
             const LinearAlgebra::distributed::Vector<double> &)
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
  DiagonalMatrix<LinearAlgebra::distributed::Vector<double>> unit_matrix;
  unit_matrix.get_vector().reinit(30);
  unit_matrix.get_vector() = 1.0;

  LinearAlgebra::distributed::Vector<double> matrix_entries(unit_matrix.m());
  for (unsigned int i = 0; i < unit_matrix.m(); ++i)
    matrix_entries(i) = i + 1;
  MyDiagonalMatrix matrix(matrix_entries);

  LinearAlgebra::distributed::Vector<double> rhs(unit_matrix.m()),
    sol(unit_matrix.m());
  rhs = 1.;

  deallog << "Solve with PreconditionIdentity: " << std::endl;
  SolverControl                                        control(30, 1e-4);
  SolverCG<LinearAlgebra::distributed::Vector<double>> solver(control);
  solver.connect(&monitor_norm);
  solver.solve(matrix, sol, rhs, PreconditionIdentity());

  deallog << "Solve with diagonal preconditioner: " << std::endl;
  sol = 0;
  solver.solve(matrix, sol, rhs, unit_matrix);
}
