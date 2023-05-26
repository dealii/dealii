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
 * Find the solution of f(x)=sin(x)=0, which is rather simple to solve, starting
 * at x=1.
 *
 * This test lets the residual function throw an exception if the step
 * length was too large from the previous point of evaluation. Let's
 * see whether PETSc SNES manages to recover by reducing the step
 * length in case of a recoverable error.
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

    const double starting_x         = 1;
    double       last_residual_eval = starting_x;

    solver.residual = [&](const VectorType &X, VectorType &F) -> void {
      deallog << "Evaluating the residual at x=" << X(0) << std::endl;

      if (std::abs(X(0) - last_residual_eval) > 1)
        {
          deallog << "Too large a step. Throwing a recoverable exception."
                  << std::endl;
          throw RecoverableUserCallbackError();
        }

      F(0) = std::sin(X(0));
      F.compress(VectorOperation::insert);

      last_residual_eval = X(0);
    };

    for (int setjac = 0; setjac < 2; setjac++)
      {
        if (setjac)
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
        x(0) = starting_x;
        x.compress(VectorOperation::insert);

        deallog << "Solving with setjac " << setjac << std::endl;
        try
          {
            solver.reinit(data);
            const auto nit = solver.solve(x);

            deallog << "Found the solution x=" << x(0) << " after " << nit
                    << " iterations." << std::endl;
          }
        catch (std::exception &exc)
          {
            deallog << exc.what() << std::endl;
          }
      }
  }
}
