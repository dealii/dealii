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

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_snes.h>

#include <cmath>

#include "../tests.h"

/**
 * Solves the nonlinear system of equations
 *
 * (x - y^3 + 1)^3 - y^3 = 0
 * x + 2y - 3 = 0
 *
 * using the PETScWrappers::NonlinearSolver class
 * that interfaces PETSc SNES solver object.
 *
 * This code tests the block interface.
 * See petsc_snes_00.cc for additional information on the callbacks
 */
using VectorType = PETScWrappers::MPI::BlockVector;
using MatrixType = PETScWrappers::MPI::BlockSparseMatrix;
using Solver     = PETScWrappers::NonlinearSolver<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  std::ofstream out("output");

  {
    out << "# Test PETSc interface (MFFD)" << std::endl;

    Solver solver;

    solver.residual = [&](const VectorType &X, VectorType &F) -> void {
      auto x = X.block(0)[0];
      auto y = X.block(1)[0];

      F.block(0)[0] = std::pow(x - std::pow(y, 3) + 1, 3) - std::pow(y, 3);
      F.block(1)[0] = x + 2 * y - 3;
      F.compress(VectorOperation::insert);
    };

    VectorType x(2, MPI_COMM_SELF, 1, 1);
    x = 0.0;

    auto nit = solver.solve(x);

    out << "#   Solution " << x[0] << ", " << x[1] << std::endl;
    out << "#   Iterations " << nit << std::endl;
  }

  {
    out << "# Test PETSc interface" << std::endl;

    Solver solver;

    solver.residual = [&](const VectorType &X, VectorType &F) -> void {
      auto x = X.block(0)[0];
      auto y = X.block(1)[0];

      F.block(0)[0] = std::pow(x - std::pow(y, 3) + 1, 3) - std::pow(y, 3);
      F.block(1)[0] = x + 2 * y - 3;
      F.compress(VectorOperation::insert);
    };

    solver.jacobian =
      [&](const VectorType &X, MatrixType &A, MatrixType &P) -> void {
      auto x    = X.block(0)[0];
      auto y    = X.block(1)[0];
      auto f0_x = 3 * std::pow(x - std::pow(y, 3) + 1, 2);
      auto f0_y = -9 * std::pow(x - std::pow(y, 3) + 1, 2) * std::pow(y, 2) -
                  3 * std::pow(y, 2);
      P.block(0, 0).set(0, 0, f0_x);
      P.block(0, 1).set(0, 0, f0_y);
      P.block(1, 0).set(0, 0, 1);
      P.block(1, 1).set(0, 0, 2);
      P.compress(VectorOperation::insert);
    };

    VectorType x(2, MPI_COMM_SELF, 1, 1);
    x = 0.0;

    // We need to create a representative matrix upfront
    MatrixType J;
    J.reinit(2, 2);
    DynamicSparsityPattern csp(1, 1);
    csp.add(0, 0);
    IndexSet indices(1);
    indices.add_range(0, 1);
    for (unsigned int row = 0; row < 2; ++row)
      for (unsigned int col = 0; col < 2; ++col)
        J.block(row, col).reinit(indices, indices, csp, MPI_COMM_SELF);
    J.collect_sizes();

    // Solve the nonlinear system and specify the matrix instance to
    // use. The matrix will be returned back to the user for assembly.
    // Note that this is currently done by wrapping PETSc's Mat objects
    // within light-weight on-the-fly constructors.
    // Any additional information stored in non standard members of J
    // will be lost, unless it can be deducted from the Mat object.
    auto nit = solver.solve(x, J);

    out << "#   Solution " << x[0] << ", " << x[1] << std::endl;
    out << "#   Iterations " << nit << std::endl;
  }

  {
    out << "# Test user interface" << std::endl;

    Solver solver;

    solver.residual = [&](const VectorType &X, VectorType &F) -> void {
      auto x = X.block(0)[0];
      auto y = X.block(1)[0];

      F.block(0)[0] = std::pow(x - std::pow(y, 3) + 1, 3) - std::pow(y, 3);
      F.block(1)[0] = x + 2 * y - 3;
      F.compress(VectorOperation::insert);
    };

    FullMatrix<double> Jinv(2, 2);

    solver.setup_jacobian = [&](const VectorType &X) -> void {
      auto x    = X.block(0)[0];
      auto y    = X.block(1)[0];
      auto f0_x = 3 * std::pow(x - std::pow(y, 3) + 1, 2);
      auto f0_y = -9 * std::pow(x - std::pow(y, 3) + 1, 2) * std::pow(y, 2) -
                  3 * std::pow(y, 2);
      FullMatrix<double> J(2, 2);
      J(0, 0) = f0_x;
      J(0, 1) = f0_y;
      J(1, 0) = 1;
      J(1, 1) = 2;
      Jinv.invert(J);
    };

    solver.solve_with_jacobian = [&](const VectorType &src,
                                     VectorType       &dst) -> void {
      dst.block(0)[0] =
        Jinv(0, 0) * src.block(0)[0] + Jinv(0, 1) * src.block(1)[0];
      dst.block(1)[0] =
        Jinv(1, 0) * src.block(0)[0] + Jinv(1, 1) * src.block(1)[0];
      dst.compress(VectorOperation::insert);
    };

    VectorType x(2, MPI_COMM_SELF, 1, 1);
    x = 0.0;

    auto nit = solver.solve(x);

    out << "#   Solution " << x[0] << ", " << x[1] << std::endl;
    out << "#   Iterations " << nit << std::endl;
  }
}
