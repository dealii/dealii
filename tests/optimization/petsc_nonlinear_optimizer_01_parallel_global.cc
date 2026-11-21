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

// On four MPI processes, solve min 1/2 sum_{i=0}^7 (x_i-(i+1))^2 +
// 1/2 (sum_{i=0}^7 x_i-36)^2. The second term couples all variables, and the
// minimizer is x_i=i+1.

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

  const MPI_Comm comm  = MPI_COMM_WORLD;
  const PetscInt first = 2 * Utilities::MPI::this_mpi_process(comm);
  VectorType     solution_vector(comm, 8, 2);
  MatrixType     hessian_matrix;
  Mat            petsc_hessian;
  PetscCallAbort(
    comm,
    MatCreateAIJ(comm, 2, 2, 8, 8, 2, nullptr, 6, nullptr, &petsc_hessian));
  static_cast<PETScWrappers::MatrixBase &>(hessian_matrix)
    .reinit(petsc_hessian);
  PetscCallAbort(comm, MatDestroy(&petsc_hessian));

  Optimizer optimizer;
  optimizer.objective = [first, comm](const VectorType &x) -> PetscReal {
    const PetscReal local_sum = x(first) + x(first + 1);
    const PetscReal sum       = Utilities::MPI::sum(local_sum, comm);
    return Utilities::MPI::sum(.5 * ((x(first) - (first + 1.)) *
                                       (x(first) - (first + 1.)) +
                                     (x(first + 1) - (first + 2.)) *
                                       (x(first + 1) - (first + 2.))),
                               comm) +
           .5 * (sum - 36.) * (sum - 36.);
  };
  optimizer.gradient = [first, comm](const VectorType &x,
                                     VectorType       &gradient_vector) {
    const PetscReal sum    = Utilities::MPI::sum(x(first) + x(first + 1), comm);
    gradient_vector(first) = x(first) - (first + 1.) + sum - 36.;
    gradient_vector(first + 1) = x(first + 1) - (first + 2.) + sum - 36.;
    gradient_vector.compress(VectorOperation::insert);
  };
  optimizer.hessian =
    [first](const VectorType &, MatrixType &hessian_matrix, MatrixType &) {
      hessian_matrix = 0.;
      for (PetscInt row = first; row < first + 2; ++row)
        for (PetscInt column = 0; column < 8; ++column)
          hessian_matrix.set(row, column, row == column ? 2. : 1.);
      hessian_matrix.compress(VectorOperation::insert);
    };

  solution_vector       = 0.;
  const auto iterations = optimizer.solve(solution_vector, hessian_matrix);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  AssertThrow(std::abs(solution_vector(first) - (first + 1.)) < 1e-4 &&
                std::abs(solution_vector(first + 1) - (first + 2.)) < 1e-4,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
