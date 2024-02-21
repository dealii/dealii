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
 * Tests exception raise for non convergence.
 */
using VectorType = PETScWrappers::MPI::Vector;
using Solver     = PETScWrappers::NonlinearSolver<VectorType>;
using real_type  = Solver::real_type;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  {
    PETScWrappers::NonlinearSolverData data;
    data.maximum_non_linear_iterations = 0;

    Solver solver(data);

    solver.residual = [&](const VectorType &X, VectorType &F) -> void {
      F.equ(2, X);
    };

    auto       commsize = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    auto       commrank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    VectorType x(MPI_COMM_WORLD, 10, commrank == commsize - 1 ? 10 : 0);
    x = 1.0;

    try
      {
        auto nit = solver.solve(x);
      }
    catch (const std::exception &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }
}
