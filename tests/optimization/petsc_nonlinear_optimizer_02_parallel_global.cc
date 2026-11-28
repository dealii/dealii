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

// On four MPI processes, solve min 1/2 sum_{i=0}^7 x_i^2, subject to the
// global equality sum_{i=0}^7 x_i=12 and x_0>=1. The equality couples all
// variables, and the minimizer is x_i=3/2.

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/optimization/petsc_nonlinear_optimizer.templates.h>

#include <cmath>

#include "../tests.h"

using VectorType = PETScWrappers::MPI::Vector;
using MatrixType = PETScWrappers::MPI::SparseMatrix;
using Optimizer  = PETScWrappers::NonlinearOptimizer<VectorType, MatrixType>;

namespace
{
  void
  reinit_constraint_matrix(MatrixType &matrix, const MPI_Comm comm)
  {
    const PetscInt local_rows =
      Utilities::MPI::this_mpi_process(comm) == 0 ? 1 : 0;
    Mat petsc_matrix;
    PetscCallAbort(
      comm,
      MatCreateAIJ(
        comm, local_rows, 2, 1, 8, 2, nullptr, 6, nullptr, &petsc_matrix));
    static_cast<PETScWrappers::MatrixBase &>(matrix).reinit(petsc_matrix);
    PetscCallAbort(comm, MatDestroy(&petsc_matrix));
  }
} // namespace



int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  const MPI_Comm comm  = MPI_COMM_WORLD;
  const PetscInt rank  = Utilities::MPI::this_mpi_process(comm);
  const PetscInt first = 2 * rank;
  VectorType     solution_vector(comm, 8, 2);
  MatrixType     hessian_matrix;
  Mat            petsc_hessian;
  PetscCallAbort(
    comm,
    MatCreateAIJ(comm, 2, 2, 8, 8, 2, nullptr, 6, nullptr, &petsc_hessian));
  static_cast<PETScWrappers::MatrixBase &>(hessian_matrix)
    .reinit(petsc_hessian);
  PetscCallAbort(comm, MatDestroy(&petsc_hessian));

  std::optional<Optimizer::Constraints> equalities;
  std::optional<Optimizer::Constraints> inequalities;
  equalities.emplace();
  equalities->values.reinit(comm, 1, rank == 0 ? 1 : 0);
  reinit_constraint_matrix(equalities->jacobian, comm);
  inequalities.emplace();
  inequalities->values.reinit(comm, 1, rank == 0 ? 1 : 0);
  reinit_constraint_matrix(inequalities->jacobian, comm);

  Optimizer optimizer;
  optimizer.objective = [first, comm](const VectorType &x) -> PetscReal {
    return Utilities::MPI::sum(.5 * (x(first) * x(first) +
                                     x(first + 1) * x(first + 1)),
                               comm);
  };
  optimizer.gradient = [first](const VectorType &x,
                               VectorType       &gradient_vector) {
    gradient_vector(first)     = x(first);
    gradient_vector(first + 1) = x(first + 1);
    gradient_vector.compress(VectorOperation::insert);
  };
  optimizer.hessian =
    [first](const VectorType &, MatrixType &hessian_matrix, MatrixType &) {
      hessian_matrix = 0.;
      hessian_matrix.set(first, first, 1.);
      hessian_matrix.set(first + 1, first + 1, 1.);
      hessian_matrix.compress(VectorOperation::insert);
    };
  optimizer.equality_constraints = [first, rank, comm](const VectorType &x,
                                                       VectorType &values) {
    const PetscReal sum = Utilities::MPI::sum(x(first) + x(first + 1), comm);
    if (rank == 0)
      values(0) = sum - 12.;
    values.compress(VectorOperation::insert);
  };
  optimizer.equality_jacobian =
    [rank](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      if (rank == 0)
        for (PetscInt column = 0; column < 8; ++column)
          matrix.set(0, column, 1.);
      matrix.compress(VectorOperation::insert);
    };
  optimizer.inequality_constraints = [rank](const VectorType &x,
                                            VectorType       &values) {
    if (rank == 0)
      values(0) = x(0) - 1.;
    values.compress(VectorOperation::insert);
  };
  optimizer.inequality_jacobian =
    [rank](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      if (rank == 0)
        matrix.set(0, 0, 1.);
      matrix.compress(VectorOperation::insert);
    };

  solution_vector = 0.;
  const auto iterations =
    optimizer.solve(solution_vector, hessian_matrix, equalities, inequalities);
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;
  AssertThrow(std::abs(solution_vector(first) - 1.5) < 1e-5 &&
                std::abs(solution_vector(first + 1) - 1.5) < 1e-5,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
