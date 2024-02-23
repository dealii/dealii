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
 * Find the solution of f(x)=sin(x)=0, which is rather simple to solve, starting
 * at x=1.
 */
using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Solver     = PETScWrappers::NonlinearSolver<VectorType, MatrixType>;
using real_type  = Solver::real_type;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  {
    PETScWrappers::NonlinearSolverData data;
    Solver                             solver(data);

    solver.residual = [&](const VectorType &X, VectorType &F) -> void {
      deallog << "Evaluating the residual at x=" << X(0) << std::endl;

      F(0) = std::sin(X(0));
      F.compress(VectorOperation::insert);
    };

    solver.jacobian =
      [&](const VectorType &X, MatrixType &A, MatrixType &P) -> void {
      deallog << "Evaluating the Jacobian at x=" << X(0) << std::endl;

      P.set(0, 0, std::cos(X(0)));
      P.compress(VectorOperation::insert);
      A.set(0, 0, std::cos(X(0)));
      A.compress(VectorOperation::insert);
    };

    auto       commsize = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    auto       commrank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    VectorType x(MPI_COMM_WORLD, 1, commrank == commsize - 1 ? 1 : 0);
    x(0) = 1.0;
    x.compress(VectorOperation::insert);

    try
      {
        const auto nit = solver.solve(x);

        deallog << "Found the solution x=" << x(0) << " after " << nit
                << " iterations." << std::endl;
      }
    catch (const std::exception &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }
}
