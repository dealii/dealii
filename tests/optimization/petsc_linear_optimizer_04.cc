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

// Solve the unconstrained linear program min 0. Every point is a minimizer, so
// the optimizer retains the initial guess (1,-2).

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_linear_optimizer.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;

int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType solution_vector(MPI_COMM_SELF, 2, 2);
  VectorType objective(MPI_COMM_SELF, 2, 2);
  solution_vector(0) = 1.;
  solution_vector(1) = -2.;
  solution_vector.compress(VectorOperation::insert);
  objective = 0.;

  PETScWrappers::LinearOptimizer<VectorType, MatrixType> optimizer;
  const auto iterations = optimizer.solve(solution_vector, objective);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  AssertThrow(std::abs(solution_vector(0) - 1.) < 1e-6 &&
                std::abs(solution_vector(1) + 2.) < 1e-6,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
