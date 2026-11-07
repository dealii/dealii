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

// Solve min 1/2 (x_0^2 + x_1^2), subject to x_0 + x_1 = 3 and x_0 >= 1.
// The unique minimizer is (3/2, 3/2).

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
  reinit_matrix(MatrixType &matrix, const PetscInt rows, const PetscInt columns)
  {
    Mat petsc_matrix;
    PetscCallAbort(
      MPI_COMM_SELF,
      MatCreateSeqAIJ(
        MPI_COMM_SELF, rows, columns, columns, nullptr, &petsc_matrix));
    static_cast<PETScWrappers::MatrixBase &>(matrix).reinit(petsc_matrix);
    PetscCallAbort(MPI_COMM_SELF, MatDestroy(&petsc_matrix));
  }
} // namespace



int
main(int argc, char **argv)
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  VectorType                            solution_vector(MPI_COMM_SELF, 2, 2);
  MatrixType                            hessian_matrix;
  Optimizer                             optimizer;
  std::optional<Optimizer::Constraints> equalities;
  std::optional<Optimizer::Constraints> inequalities;

  reinit_matrix(hessian_matrix, 2, 2);
  equalities.emplace();
  equalities->values.reinit(MPI_COMM_SELF, 1, 1);
  reinit_matrix(equalities->jacobian, 1, 2);
  inequalities.emplace();
  inequalities->values.reinit(MPI_COMM_SELF, 1, 1);
  reinit_matrix(inequalities->jacobian, 1, 2);

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
  optimizer.equality_constraints = [](const VectorType &x, VectorType &values) {
    values(0) = x(0) + x(1) - 3.;
    values.compress(VectorOperation::insert);
  };
  optimizer.equality_jacobian =
    [](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      matrix.set(0, 0, 1.);
      matrix.set(0, 1, 1.);
      matrix.compress(VectorOperation::insert);
    };
  optimizer.inequality_constraints = [](const VectorType &x,
                                        VectorType       &values) {
    values(0) = x(0) - 1.;
    values.compress(VectorOperation::insert);
  };
  optimizer.inequality_jacobian =
    [](const VectorType &, MatrixType &matrix, MatrixType &) {
      matrix = 0.;
      matrix.set(0, 0, 1.);
      matrix.compress(VectorOperation::insert);
    };

  solution_vector = 0.;
  const auto iterations =
    optimizer.solve(solution_vector, hessian_matrix, equalities, inequalities);
  deallog << "The optimizer found the solution [" << solution_vector(0) << ','
          << solution_vector(1) << "]." << std::endl;
  deallog << "The optimizer found the solution in " << iterations
          << " iterations." << std::endl;

  // Now double-check that that is actually correct:
  AssertThrow(std::abs(solution_vector(0) - 1.5) < 1e-5 &&
                std::abs(solution_vector(1) - 1.5) < 1e-5,
              ExcInternalError());
  deallog << "OK" << std::endl;
}
