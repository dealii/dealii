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
 *
 * This test lets the jacobian callbacks throw an exception if the step
 * length was too large from the previous point of evaluation. SNES is not able
 * able to recover from such errors.
 * Here we test that we propagate the proper converged reason up to the solve
 * callpoint in case of a recoverable exception.
 * We also test unrecoverable errors handling.
 */
using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Solver     = PETScWrappers::NonlinearSolver<VectorType, MatrixType>;
using real_type  = Solver::real_type;

void
run_test(int testcase, bool recoverable)
{
  PETScWrappers::NonlinearSolverData data;
  Solver                             solver(data);

  const double starting_x         = 1;
  double       last_residual_eval = starting_x;
  bool         throw_exception    = false;
  bool         throw_recoverable  = recoverable;

  solver.residual = [&](const VectorType &X, VectorType &F) -> void {
    deallog << "Evaluating the residual at x=" << X(0) << std::endl;

    F(0) = std::sin(X(0));
    F.compress(VectorOperation::insert);
    if (std::abs(X(0) - last_residual_eval) > 1)
      throw_exception = true;

    last_residual_eval = X(0);
  };

  double invJ;

  if (testcase)
    {
      solver.setup_jacobian = [&](const VectorType &X) -> void {
        deallog << "Setup Jacobian at x=" << X(0) << std::endl;

        invJ = 1. / std::cos(X(0));
        if (throw_exception && testcase == 1)
          {
            if (throw_recoverable)
              {
                deallog
                  << "Throwing a recoverable exception from setup_jacobian."
                  << std::endl;
                throw RecoverableUserCallbackError();
              }
            else
              {
                deallog
                  << "Throwing a unrecoverable exception from setup_jacobian."
                  << std::endl;
                throw std::runtime_error(
                  "Unrecoverable error in setup_jacobian");
              }
          }
      };

      solver.solve_with_jacobian = [&](const VectorType &rhs,
                                       VectorType       &sol) -> void {
        deallog << "Solve Jacobian" << std::endl;

        if (throw_exception && testcase == 2)
          {
            if (throw_recoverable)
              {
                deallog
                  << "Throwing a recoverable exception from solve_with_jacobian."
                  << std::endl;
                throw RecoverableUserCallbackError();
              }
            else
              {
                deallog
                  << "Throwing a unrecoverable exception from solve_with_jacobian."
                  << std::endl;
                throw std::runtime_error(
                  "Unrecoverable error in solve_with_jacobian");
              }
          }
        else
          {
            sol(0) = invJ * rhs(0);
            sol.compress(VectorOperation::insert);
          }
      };
    }
  else
    {
      solver.jacobian =
        [&](const VectorType &X, MatrixType &A, MatrixType &P) -> void {
        deallog << "Evaluating the Jacobian at x=" << X(0) << std::endl;

        if (throw_exception)
          {
            if (throw_recoverable)
              {
                deallog << "Throwing a recoverable exception from jacobian."
                        << std::endl;
                throw RecoverableUserCallbackError();
              }
            else
              {
                deallog << "Throwing a unrecoverable exception from jacobian."
                        << std::endl;
                throw std::runtime_error("Unrecoverable error in jacobian");
              }
          }
        else
          {
            P.set(0, 0, std::cos(X(0)));
            P.compress(VectorOperation::insert);
          }
      };
    }

  auto       commsize = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  auto       commrank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  VectorType x(MPI_COMM_WORLD, 1, commrank == commsize - 1 ? 1 : 0);
  x(0) = starting_x;
  x.compress(VectorOperation::insert);

  // This test triggers false positives in FPE trapping for some versions of
  // PETSc
#if DEAL_II_PETSC_VERSION_LT(3, 19, 2) && defined(DEBUG) && \
  defined(DEAL_II_HAVE_FP_EXCEPTIONS)
  PetscErrorCode ierr = PetscFPTrapPush(PETSC_FP_TRAP_OFF);
  (void)ierr;
#endif

  deallog << "Running testcase " << testcase << std::endl;
  try
    {
      const auto nit = solver.solve(x);

      deallog << "Found the solution x=" << x(0) << " after " << nit
              << " iterations." << std::endl;
    }
  catch (const std::exception &exc)
    {
      deallog << "Caught exception." << std::endl;
      deallog << exc.what() << std::endl;
    }
}

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  for (int i = 0; i < 3; i++)
    run_test(i, true); // recoverable errors
  for (int i = 0; i < 3; i++)
    run_test(i, false); // unrecoverable errors
}
