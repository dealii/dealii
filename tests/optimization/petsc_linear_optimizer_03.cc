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

// Solve min x_0+2x_1, subject to x_0>=1 and x_1>=2. The minimizer is (1,2).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Bounds =
  PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>::Bounds;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType solution_vector(MPI_COMM_SELF, 2, 2);
  VectorType objective(MPI_COMM_SELF, 2, 2);
  Bounds     bounds;
  bounds.lower.reinit(MPI_COMM_SELF, 2, 2);
  bounds.upper.reinit(MPI_COMM_SELF, 2, 2);

  solution_vector(0) = 2.;
  solution_vector(1) = 3.;
  solution_vector.compress(VectorOperation::insert);
  objective(0) = 1.;
  objective(1) = 2.;
  objective.compress(VectorOperation::insert);
  bounds.lower(0) = 1.;
  bounds.lower(1) = 2.;
  bounds.upper(0) = PETSC_INFINITY;
  bounds.upper(1) = PETSC_INFINITY;
  bounds.lower.compress(VectorOperation::insert);
  bounds.upper.compress(VectorOperation::insert);

  PETScWrappers::LinearOptimizer<VectorType, MatrixType> optimizer;
  const auto iterations = optimizer.solve(
    solution_vector, objective, nullptr, nullptr, nullptr, nullptr, bounds);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  AssertThrow(std::abs(solution_vector(0) - 1.) < 1e-6 &&
                std::abs(solution_vector(1) - 2.) < 1e-6,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
