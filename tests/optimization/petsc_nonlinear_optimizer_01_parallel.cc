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

// On four MPI processes, solve the eight-variable unconstrained problem
// min 1/2 sum_{r=0}^3 ((x_{2r}-1)^2+(x_{2r+1}-2)^2). This is four
// independent quadratic pairs, with minimizer x_{2r}=1 and x_{2r+1}=2.

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

  VectorType solution_vector(comm, 8, 2);
  MatrixType hessian_matrix;
  Mat        petsc_hessian;
  PetscCallAbort(
    comm,
    MatCreateAIJ(comm, 2, 2, 8, 8, 1, nullptr, 0, nullptr, &petsc_hessian));
  static_cast<PETScWrappers::MatrixBase &>(hessian_matrix)
    .reinit(petsc_hessian);
  PetscCallAbort(comm, MatDestroy(&petsc_hessian));

  Optimizer optimizer;
  optimizer.objective = [first, comm](const VectorType &x) -> PetscReal {
    const PetscReal local_value =
      .5 * ((x(first) - 1.) * (x(first) - 1.) +
            (x(first + 1) - 2.) * (x(first + 1) - 2.));
    return Utilities::MPI::sum(local_value, comm);
  };
  optimizer.gradient = [first](const VectorType &x,
                               VectorType       &gradient_vector) {
    gradient_vector(first)     = x(first) - 1.;
    gradient_vector(first + 1) = x(first + 1) - 2.;
    gradient_vector.compress(VectorOperation::insert);
  };
  optimizer.hessian =
    [first](const VectorType &, MatrixType &hessian_matrix, MatrixType &) {
      hessian_matrix = 0.;
      hessian_matrix.set(first, first, 1.);
      hessian_matrix.set(first + 1, first + 1, 1.);
      hessian_matrix.compress(VectorOperation::insert);
    };

  solution_vector       = 0.;
  const auto iterations = optimizer.solve(solution_vector, hessian_matrix);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  AssertThrow(std::abs(solution_vector(first) - 1.) < 1e-6 &&
                std::abs(solution_vector(first + 1) - 2.) < 1e-6,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
