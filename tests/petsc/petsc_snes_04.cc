//-----------------------------------------------------------
//
//    Copyright (C) 2023 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

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
    catch (std::exception &exc)
      {
        deallog << exc.what() << std::endl;
      }
  }
}
