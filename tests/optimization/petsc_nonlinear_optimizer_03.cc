// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Solve min 1/2 (x_0^2 + x_1^2), subject to x_0 >= 1. The optimizer rejects
// the infeasible initial guess (0, 0) and converges from (1, 1) to (1, 0).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType solution_vector(MPI_COMM_SELF, 2, 2);
  MatrixType hessian_matrix;
  Mat        petsc_hessian;
  PetscCallAbort(
    MPI_COMM_SELF,
    MatCreateSeqAIJ(MPI_COMM_SELF, 2, 2, 2, nullptr, &petsc_hessian));
  static_cast<PETScWrappers::MatrixBase &>(hessian_matrix)
    .reinit(petsc_hessian);
  PetscCallAbort(MPI_COMM_SELF, MatDestroy(&petsc_hessian));

  std::optional<Optimizer::Bounds> bounds;
  bounds.emplace();
  bounds->lower.reinit(MPI_COMM_SELF, 2, 2);
  bounds->upper.reinit(MPI_COMM_SELF, 2, 2);
  bounds->lower(0) = 1.;
  bounds->lower(1) = PETSC_NINFINITY;
  bounds->upper(0) = PETSC_INFINITY;
  bounds->upper(1) = PETSC_INFINITY;
  bounds->lower.compress(VectorOperation::insert);
  bounds->upper.compress(VectorOperation::insert);

  Optimizer optimizer;
  optimizer.objective = [](const VectorType &x) -> PetscReal {
    return .5 * (x(0) * x(0) + x(1) * x(1));
  };
  optimizer.gradient = [](const VectorType &x, VectorType &gradient_vector) {
    gradient_vector(0) = x(0);
    gradient_vector(1) = x(1);
    gradient_vector.compress(VectorOperation::insert);
  };
  optimizer.hessian =
    [](const VectorType &, MatrixType &hessian_matrix, MatrixType &) {
      hessian_matrix = 0.;
      hessian_matrix.set(0, 0, 1.);
      hessian_matrix.set(1, 1, 1.);
      hessian_matrix.compress(VectorOperation::insert);
    };

  solution_vector                        = 0.;
  bool infeasible_initial_guess_rejected = false;
  try
    {
      optimizer.solve(
        solution_vector, hessian_matrix, std::nullopt, std::nullopt, bounds);
    }
  catch (const ExceptionBase &)
    {
      infeasible_initial_guess_rejected = true;
    }
  AssertThrow(infeasible_initial_guess_rejected, ExcInternalError());

  solution_vector(0) = 1.;
  solution_vector(1) = 1.;
  solution_vector.compress(VectorOperation::insert);
  const auto iterations = optimizer.solve(
    solution_vector, hessian_matrix, std::nullopt, std::nullopt, bounds);
  deallog << "The optimizer found the solution [" << solution_vector(0) << ','
          << solution_vector(1) << "]." << std::endl;
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  // Now double-check that that is actually correct:
  AssertThrow(std::abs(solution_vector(0) - 1.) < 1e-5 &&
                std::abs(solution_vector(1)) < 1e-5,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
