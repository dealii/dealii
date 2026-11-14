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

// On four MPI processes, solve min 1/2 sum_{i=0}^7 x_i^2, subject to
// x_{2r}+x_{2r+1}=3 and x_{2r}>=1 for r=0,...,3. This is four independent
// constrained pairs, each with minimizer (3/2,3/2).

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
    Mat petsc_matrix;
    PetscCallAbort(
      comm,
      MatCreateAIJ(comm, 1, 2, 4, 8, 2, nullptr, 0, nullptr, &petsc_matrix));
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
    MatCreateAIJ(comm, 2, 2, 8, 8, 1, nullptr, 0, nullptr, &petsc_hessian));
  static_cast<PETScWrappers::MatrixBase &>(hessian_matrix)
    .reinit(petsc_hessian);
  PetscCallAbort(comm, MatDestroy(&petsc_hessian));

  std::optional<Optimizer::Constraints> equalities;
  std::optional<Optimizer::Constraints> inequalities;
  equalities.emplace();
  equalities->values.reinit(comm, 4, 1);
  reinit_constraint_matrix(equalities->jacobian, comm);
  inequalities.emplace();
  inequalities->values.reinit(comm, 4, 1);
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
  optimizer.equality_constraints = [first, rank](const VectorType &x,
                                                 VectorType       &values) {
    values(rank) = x(first) + x(first + 1) - 3.;
    values.compress(VectorOperation::insert);
  };
  optimizer.equality_jacobian =
    [first, rank](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      matrix.set(rank, first, 1.);
      matrix.set(rank, first + 1, 1.);
      matrix.compress(VectorOperation::insert);
    };
  optimizer.inequality_constraints = [first, rank](const VectorType &x,
                                                   VectorType       &values) {
    values(rank) = x(first) - 1.;
    values.compress(VectorOperation::insert);
  };
  optimizer.inequality_jacobian =
    [first, rank](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      matrix.set(rank, first, 1.);
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
