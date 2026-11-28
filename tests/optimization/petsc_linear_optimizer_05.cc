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

// Solve min sum_{i=0}^{19} x_i, subject to x_i>=0 and x_{i+1}-x_i>=1 for
// i=0,...,18. The unique minimizer is x_i=i. The initial guess is
// x_i=i+10+U(0,100).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.h>

#include <cmath>
#include <random>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  constexpr PetscInt n_variables    = 20;
  constexpr PetscInt n_inequalities = n_variables - 1;
  const MPI_Comm     comm           = MPI_COMM_SELF;

  VectorType solution_vector(comm, n_variables, n_variables);
  VectorType objective(comm, n_variables, n_variables);
  VectorType inequality_right_hand_side(comm, n_inequalities, n_inequalities);
  MatrixType inequality_matrix;
  Mat        petsc_matrix;
  PetscCallAbort(comm,
                 MatCreateAIJ(comm,
                              n_inequalities,
                              n_variables,
                              n_inequalities,
                              n_variables,
                              2,
                              nullptr,
                              0,
                              nullptr,
                              &petsc_matrix));
  static_cast<PETScWrappers::MatrixBase &>(inequality_matrix)
    .reinit(petsc_matrix);
  PetscCallAbort(comm, MatDestroy(&petsc_matrix));

  for (PetscInt index = 0; index < n_variables; ++index)
    {
      std::mt19937                           generator(index);
      std::uniform_real_distribution<double> distribution(0., 100.);
      solution_vector(index) = index + 10. + distribution(generator);
      objective(index)       = 1.;
    }
  solution_vector.compress(VectorOperation::insert);
  objective.compress(VectorOperation::insert);

  for (PetscInt row = 0; row < n_inequalities; ++row)
    {
      inequality_matrix.set(row, row, -1.);
      inequality_matrix.set(row, row + 1, 1.);
      inequality_right_hand_side(row) = 1.;
    }
  inequality_matrix.compress(VectorOperation::insert);
  inequality_right_hand_side.compress(VectorOperation::insert);

  std::optional<Optimizer::Bounds> bounds;
  bounds.emplace();
  bounds->lower.reinit(comm, n_variables, n_variables);
  bounds->upper.reinit(comm, n_variables, n_variables);
  for (PetscInt index = 0; index < n_variables; ++index)
    {
      bounds->lower(index) = 0.;
      bounds->upper(index) = PETSC_INFINITY;
    }
  bounds->lower.compress(VectorOperation::insert);
  bounds->upper.compress(VectorOperation::insert);

  PETScWrappers::LinearOptimizer<VectorType, MatrixType> solver;
  const auto iterations = solver.solve(solution_vector,
                                       objective,
                                       nullptr,
                                       nullptr,
                                       &inequality_matrix,
                                       &inequality_right_hand_side,
                                       bounds);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  for (PetscInt index = 0; index < n_variables; ++index)
    AssertThrow(std::abs(solution_vector(index) - index) < 1e-6,
                ExcInternalError());
  deallog << "OK" << std::endl;
}
