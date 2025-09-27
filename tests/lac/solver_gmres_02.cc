// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Check the path of SolverGMRES with dealii::Vector and
// LinearAlgebra::OrthogonalizationStrategy::classical_gram_schmidt. Apart
// from the vector type, this test case is the same as solver_gmres_01.


#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


struct MyDiagonalMatrix
{
  MyDiagonalMatrix(const Vector<double> &diagonal)
    : diagonal(diagonal)
  {}

  void
  vmult(Vector<double> &dst, const Vector<double> &src) const
  {
    dst = src;
    dst.scale(diagonal);
  }

  void
  vmult(BlockVector<double> &dst, const BlockVector<double> &src) const
  {
    dst = src;
    for (unsigned int b = 0; b < src.n_blocks(); ++b)
      dst.block(0).scale(diagonal);
  }

  unsigned int
  m() const
  {
    return diagonal.size();
  }

  const Vector<double> &diagonal;
};



SolverControl::State
monitor_norm(const unsigned int iteration,
             const double       check_value,
             const Vector<double> &)
{
  deallog << "   estimated residual at iteration " << iteration << ": "
          << check_value << std::endl;
  return SolverControl::success;
}



SolverControl::State
monitor_norm_block(const unsigned int iteration,
                   const double       check_value,
                   const BlockVector<double> &)
{
  deallog << "   estimated residual at iteration " << iteration << ": "
          << check_value << std::endl;
  return SolverControl::success;
}


int
main()
{
  initlog();

  // Create diagonal matrix with entries between 1 and 30
  Vector<double> matrix_entries_(30);
  matrix_entries_ = 1.0;
  MyDiagonalMatrix unit_matrix(matrix_entries_);

  Vector<double> matrix_entries(unit_matrix.m());
  for (unsigned int i = 0; i < unit_matrix.m(); ++i)
    matrix_entries(i) = i + 1;
  MyDiagonalMatrix matrix(matrix_entries);

  Vector<double> rhs(unit_matrix.m());
  Vector<double> sol(unit_matrix.m());
  rhs = 1.;

  {
    deallog << "Solve with PreconditionIdentity: " << std::endl;
    SolverControl                               control(40, 1e-4);
    SolverGMRES<Vector<double>>::AdditionalData data3(6);
    data3.orthogonalization_strategy =
      LinearAlgebra::OrthogonalizationStrategy::classical_gram_schmidt;
    SolverGMRES<Vector<double>> solver(control, data3);
    solver.connect(&monitor_norm);
    solver.solve(matrix, sol, rhs, PreconditionIdentity());

    deallog << "Solve with diagonal preconditioner: " << std::endl;
    sol = 0;
    solver.solve(matrix, sol, rhs, unit_matrix);
  }

  {
    BlockVector<double> rhs_(1);
    BlockVector<double> sol_(1);
    rhs_.block(0) = rhs;
    sol_.block(0) = sol;

    deallog << "Solve with PreconditionIdentity: " << std::endl;
    SolverControl                                    control(40, 1e-4);
    SolverGMRES<BlockVector<double>>::AdditionalData data3(6);
    data3.orthogonalization_strategy =
      LinearAlgebra::OrthogonalizationStrategy::classical_gram_schmidt;
    SolverGMRES<BlockVector<double>> solver(control, data3);
    solver.connect(&monitor_norm_block);
    solver.solve(matrix, sol_, rhs_, PreconditionIdentity());

    deallog << "Solve with diagonal preconditioner: " << std::endl;
    sol_ = 0;
    solver.solve(matrix, sol_, rhs_, unit_matrix);
  }
}
